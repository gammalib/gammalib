#! /usr/bin/env python
# ==========================================================================
# This script tests the elliptical Npred computation. The script requires
# matplotlib to be installed.
#
# --------------------------------------------------------------------------
#
# Copyright (C) 2014 Juergen Knoedlseder
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ==========================================================================
import gammalib
import math


# ================= #
# Compute arclength #
# ================= #
def cta_roi_arclength(rad, dist, roi):
    """
    Returns length of circular arc within circular ROI

    Parameters:
     rad  - Circle radius in radians (<pi)
     dist - Circle centre distance to ROI centre (<pi)
     roi  - Radius of ROI in radians
    """
    # Initialise
    arclength = 0.0

    # Handle special case that circle centre matches ROI centre
    if dist == 0.0:
        if rad > roi:
            arclength = 0.0
        else:
            arclength = 2.0 * math.pi

    # ... otherwise circle and ROI centres are not identical
    else:

        # Handle special case that we evaluate exactly at the circle
        # centre. In this case we have in fact a point, and if this point
        # falls within the ROI it has a formal arclength of 2pi.
        if rad == 0.0:
            if dist > roi:
                arclength = 0.0
            else:
                arclength = 2.0 * math.pi

        # ... otherwise we have to handle the general case
        else:
            d = roi - dist
            if -rad >= d:
                arclength = 0.0
            elif rad <= d:
                arclength = 2.0 * math.pi
            else:
                cosrad  = math.cos(rad)
                sinrad  = math.sin(rad)
                cosroi  = math.cos(roi)
                cosdist = math.cos(dist)
                sindist = math.sin(dist)
                cosang  = (cosroi - cosdist * cosrad) / (sindist * sinrad)
                arclength = 2.0 * gammalib.acos(cosang)

    # Return
    return arclength


# ============== #
# Limit interval #
# ============== #
def limit_omega(min, max, domega, log=True):
    """
    """
    # Compute only if domega < pi
    if (domega < gammalib.pi):

        # Save initial interval
        min_initial = min
        max_initial = max

        # Initialise intervals
        intervals = []
        
        # Limit each of the intervals to [omega_min, omega_max]
        omega_min       = -domega
        omega_max       = +domega
        omega_min_plus  = omega_min + gammalib.twopi
        omega_max_plus  = omega_max + gammalib.twopi
        omega_min_minus = omega_min - gammalib.twopi
        omega_max_minus = omega_max - gammalib.twopi

        # Check for overlaps
        in_unshifted = False
        in_plus      = False
        in_minus     = False
        if (max_initial > omega_min and min_initial < omega_max):
            min = min_initial
            max = max_initial
            in_unshifted = True
            if (min < omega_min):
                min = omega_min
            if (max > omega_max):
                max = omega_max
            intervals.append({'min': min, 'max': max})
        if (max_initial > omega_min_plus and min_initial < omega_max_plus):
            min = min_initial
            max = max_initial
            in_plus = True
            if (min < omega_min_plus):
                min = omega_min_plus
            if (max > omega_max_plus):
                max = omega_max_plus
            intervals.append({'min': min, 'max': max})
        if (max_initial > omega_min_minus and min_initial < omega_max_minus):
            min = min_initial
            max = max_initial
            in_minus = True
            if (min < omega_min_minus):
                min = omega_min_minus
            if (max > omega_max_minus):
                max = omega_max_minus
            intervals.append({'min': min, 'max': max})

        # Debug
        if log:
            print("[%7.2f,%7.2f] => [%7.2f,%7.2f] (%i[%7.2f,%7.2f]  %i[%7.2f,%7.2f]  %i[%7.2f,%7.2f])" %
              (min_initial*gammalib.rad2deg, max_initial*gammalib.rad2deg, \
               min*gammalib.rad2deg, max*gammalib.rad2deg, \
               in_unshifted, \
               omega_min*gammalib.rad2deg, omega_max*gammalib.rad2deg, \
               in_plus, \
               omega_min_plus*gammalib.rad2deg, omega_max_plus*gammalib.rad2deg, \
               in_minus, \
               omega_min_minus*gammalib.rad2deg, omega_max_minus*gammalib.rad2deg))

    # ... just copy over the input interval
    else:
        intervals = [{'min': min, 'max': max}]

    # Return
    return intervals


# ======================= #
# Show ellptical response #
# ======================= #
def show_response(centre, semimajor, semiminor, posang, roiDir, pntDir, roi_radius, psf_max_radius, rhos, log=False):
    """
    Show ellptical response
    """
    # Create figure
    plt.figure(1)
    plt.title("Elliptical Npred computation")

    # Convert arguments to radians
    roi_radius     *= gammalib.deg2rad
    psf_max_radius *= gammalib.deg2rad
    radius_roi      = roi_radius + psf_max_radius

    # Compute some stuff
    cos_pa       = math.cos(posang * gammalib.deg2rad)
    sin_pa       = math.sin(posang * gammalib.deg2rad)
    rho_roi      = roiDir.dist(centre)
    theta_max    = semimajor * gammalib.deg2rad
    posangle_roi = centre.posang(roiDir) # omega0

    # Rho limits
    rho_min = 0.0
    rho_max = theta_max
    if (rho_roi > radius_roi):
        rho_min = rho_roi - radius_roi

    # Setup ellipse as vector
    x = []
    y = []
    for angle in range(0, 360):
        radians     = angle * gammalib.deg2rad
        cos_radians = math.cos(radians)
        sin_radians = math.sin(radians)
        x.append(semimajor * cos_radians * sin_pa - semiminor * sin_radians * cos_pa)
        y.append(semimajor * cos_radians * cos_pa + semiminor * sin_radians * sin_pa)

    # Plot ellipse
    plt.plot(x, y, 'r-')
    plt.xlim([ 2.0*semimajor,-2.0*semimajor])
    plt.ylim([-2.0*semimajor, 2.0*semimajor])
    plt.plot([0.0], [0.0], 'ro')

    # Plot ROI centre
    x_roi = roiDir.ra_deg()  - centre.ra_deg()
    y_roi = roiDir.dec_deg() - centre.dec_deg()
    plt.plot([x_roi], [y_roi], 'bo')

    # Plot ROI circle
    x = []
    y = []
    for angle in range(0, 360):
        radians     = angle * gammalib.deg2rad
        cos_radians = math.cos(radians)
        sin_radians = math.sin(radians)
        x.append(roi_radius * cos_radians * gammalib.rad2deg + x_roi)
        y.append(roi_radius * sin_radians * gammalib.rad2deg + y_roi)
    plt.plot(x, y, 'b-.')

    # Plot ROI+Psf circle
    x = []
    y = []
    for angle in range(0, 360):
        radians     = angle * gammalib.deg2rad
        cos_radians = math.cos(radians)
        sin_radians = math.sin(radians)
        x.append(radius_roi * cos_radians * gammalib.rad2deg + x_roi)
        y.append(radius_roi * sin_radians * gammalib.rad2deg + y_roi)
    plt.plot(x, y, 'b-')

    # Print some stuff
    print("Distance model-ROI (rho_roi) ............: %.2f deg" % (rho_roi*gammalib.rad2deg))
    print("Rho interval ............................: %.2f-%.2f deg" % (rho_min*gammalib.rad2deg,rho_max*gammalib.rad2deg))

    # Continue only if rho range is valid
    if (rho_max > rho_min):

        # Plot rhos
        for rho in rhos:

            # Print header
            if log:
                print("Rho ..................................: %.2f deg" % (rho))

            # Compute some stuff
            rho_rad = rho * gammalib.deg2rad

            # Compute half length of arc that lies within Psf validity circle
            domega = 0.5 * cta_roi_arclength(rho_rad, rho_roi, radius_roi)

            # Print domega
            if log:
                print("  domega .............................: %.2f deg" % (domega*gammalib.rad2deg))

            #
            if (domega > 0.0):

                # Full containment within ellipse
                if (rho < semiminor):
            
                    # Compute phi integration range
                    omega_min = -domega
                    omega_max = +domega

                    # Print
                    if log:
                        print("  Omega range ........................: %.2f - %.2f deg" % (omega_min*gammalib.rad2deg,omega_max*gammalib.rad2deg))

                    # Plot segment
                    n    = 100
                    x    = []
                    y    = []
                    step = (omega_max - omega_min)/n
                    for i in range(n+1):
                        radians = i*step + omega_min
                        cos_radians = math.sin(radians+posangle_roi)
                        sin_radians = math.cos(radians+posangle_roi)
                        x.append(rho * cos_radians)
                        y.append(rho * sin_radians)
                    plt.plot(x, y, 'k-')

                # ... otherwise compute the arcs
                else:

                    # Azimuth angle of intersection points between ellipse and
                    # circle with radius of rho
                    arg1        = 1.0 - (semiminor*semiminor) / (rho*rho)
                    arg2        = 1.0 - (semiminor*semiminor) / (semimajor*semimajor)
                    omega_width = math.acos(math.sqrt(arg1/arg2))

                    # Print Omega width
                    if log:
                        print("  omega_width ........................: %.2f deg" % (omega_width*gammalib.rad2deg))

                    # 
                    if (omega_width > 0.0):
                
                        # Compute Omega intervals. Recall that Omega is with respect
                        # to connecting line between model and observed photon. Make
                        # sure that omega_0 is within [-pi,pi] thus that the omega
                        # intervals are between [-2pi,2pi]
                        omega_0 = posang * gammalib.deg2rad - posangle_roi
                        if (omega_0 > gammalib.pi):
                            omega_0 -= gammalib.pi
                        elif (omega_0 < -gammalib.pi):
                            omega_0 += gammalib.pi

                        # Compute intervals
                        omega1_min = omega_0    - omega_width
                        omega1_max = omega_0    + omega_width
                        omega2_min = omega1_min + gammalib.pi
                        omega2_max = omega1_max + gammalib.pi

                        # Logging
                        if log:
                            print("  Omega range 1 (initial) ............: %.2f - %.2f deg" % (omega1_min*gammalib.rad2deg,omega1_max*gammalib.rad2deg))
                            print("  Omega range 2 (initial) ............: %.2f - %.2f deg" % (omega2_min*gammalib.rad2deg,omega2_max*gammalib.rad2deg))

                        # Limit intervals
                        intervals  = limit_omega(omega1_min, omega1_max, domega)
                        intervals2 = limit_omega(omega2_min, omega2_max, domega)
                        intervals.extend(intervals2)

                        # Print results
                        if log:
                            for interval in intervals:
                                print("  Omega interval .....................: %.2f - %.2f deg" % (interval['min']*gammalib.rad2deg,interval['max']*gammalib.rad2deg))

                        # Draw Omega range lines and circle segments
                        for interval in intervals:
                            omega_min = interval['min']
                            omega_max = interval['max']
                            if (omega_min < omega_max):
                                n    = 100
                                x    = []
                                y    = []
                                step = (omega_max - omega_min)/n
                                for i in range(n+1):
                                    radians = i*step + omega_min
                                    cos_radians = math.sin(radians+posangle_roi)
                                    sin_radians = math.cos(radians+posangle_roi)
                                    x.append(rho * cos_radians)
                                    y.append(rho * sin_radians)
                                plt.plot(x, y, 'k-')

    # Show plot
    plt.show()

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    Test elliptical response computation
    """
    # Dump header
    print("")
    print("*************************************")
    print("* Test elliptical Npred computation *")
    print("*************************************")

    # Check if matplotlib is available
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Matplotlib is not installed on your system. Abort.")
        sys.exit()

    # Set Parameters
    roi_radius = 0.5
    delta_max  = 0.2

    # Set model parameters
    centre = gammalib.GSkyDir()
    centre.radec_deg(0.0, 0.0)
    semimajor = 1.0
    semiminor = 0.50
    posang    = 310.0

    # Set ROI direction
    roiDir = gammalib.GSkyDir()
    #roiDir.radec_deg(-0.4, -0.3)
    roiDir.radec_deg(+0.4, -0.3)
    #roiDir.radec_deg(+0.4, +0.3)

    # Set pointing direction
    pntDir = gammalib.GSkyDir()
    pntDir.radec_deg(0.1, -1.0)

    # Set rhos
    rhos = [i*0.02 for i in range(51)]

    # Show elliptical response geometry
    show_response(centre, semimajor, semiminor, posang, roiDir, pntDir, roi_radius, delta_max, rhos)
    
