#! /usr/bin/env python
# ==========================================================================
# This script tests the elliptical response computation. The script requires
# matplotlib to be installed.
#
# --------------------------------------------------------------------------
#
# Copyright (C) 2013-2014 Juergen Knoedlseder
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


# ======================= #
# Show ellptical response #
# ======================= #
def show_response(centre, semimajor, semiminor, posang, obsDir, pntDir, rhos):
    """
    Show ellptical response
    """
    # Create figure
    plt.figure(1)
    plt.title("Elliptical response")

    # Set some constants
    delta_max = 0.5 * gammalib.deg2rad

    # Put posang in [-pi,pi]
    #if (posang > gammalib.pi):
    #    posang -= gammalib.pi
    #elif (posang < -gammalib.pi):
    #    posang += gammalib.pi

    # Compute some stuff
    cos_pa   = math.cos(posang * gammalib.deg2rad)
    sin_pa   = math.sin(posang * gammalib.deg2rad)
    zeta     = centre.dist(obsDir)
    obsOmega = centre.posang(obsDir)
    eta      = pntDir.dist(obsDir)
    lamda    = centre.dist(pntDir)
    src_max  = semimajor * gammalib.deg2rad

    # Omega_0
    omega0 = 0.0
    denom  = math.sin(lamda*gammalib.deg2rad) * math.sin(zeta*gammalib.deg2rad)
    if (denom != 0.0):
        arg    = (math.cos(eta*gammalib.deg2rad) - math.cos(lamda*gammalib.deg2rad) * \
                  math.cos(zeta*gammalib.deg2rad))/denom;
        omega0 = gammalib.acos(arg)
    pa_pnt = centre.posang(pntDir)

    # Rho limits
    rho_min = 0.0
    rho_max = zeta + delta_max
    if (zeta > delta_max):
        rho_min = zeta - delta_max
    if rho_max > src_max:
        rho_max = src_max

    # Setup ellipse as vector
    ellipse_x = []
    ellipse_y = []
    for angle in range(0, 360):
        radians     = angle * gammalib.deg2rad
        cos_radians = math.cos(radians)
        sin_radians = math.sin(radians)
        x = semimajor * cos_radians * sin_pa - semiminor * sin_radians * cos_pa
        y = semimajor * cos_radians * cos_pa + semiminor * sin_radians * sin_pa
        ellipse_x.append(x)
        ellipse_y.append(y)

    # Plot ellipse
    plt.plot(ellipse_x, ellipse_y, 'r-')
    plt.xlim([ 1.5*semimajor,-1.5*semimajor])
    plt.ylim([-1.5*semimajor, 1.5*semimajor])
    plt.plot([0.0], [0.0], 'ro')

    # Plot observed direction
    x_obs = obsDir.ra_deg()  - centre.ra_deg()
    y_obs = obsDir.dec_deg() - centre.dec_deg()
    plt.plot([x_obs], [y_obs], 'bo')

    # Draw connecting line model-photon'
    n    = 10
    step = zeta/n * gammalib.rad2deg
    x    = []
    y    = []
    for i in range(n+1):
        x.append(i*step * math.sin(obsOmega))
        y.append(i*step * math.cos(obsOmega))
    plt.plot(x, y, 'b-')

    # Plot Psf validity circle
    circle_x = []
    circle_y = []
    for angle in range(0, 360):
        radians     = angle * gammalib.deg2rad
        cos_radians = math.cos(radians)
        sin_radians = math.sin(radians)
        x = delta_max * cos_radians * gammalib.rad2deg + x_obs
        y = delta_max * sin_radians * gammalib.rad2deg + y_obs
        circle_x.append(x)
        circle_y.append(y)
    plt.plot(circle_x, circle_y, 'b-')

    # Plot pointing direction
    x = pntDir.ra_deg()  - centre.ra_deg()
    y = pntDir.dec_deg() - centre.dec_deg()
    plt.plot([x], [y], 'go')

    # Print some stuff
    print("Distance p'-model (zeta) .............: %.2f deg" % (zeta*gammalib.rad2deg))
    print("PA of p' from model (obsOmega)........: %.2f deg" % (obsOmega*gammalib.rad2deg))
    print("Distance p'-pointing (eta)............: %.2f deg" % (eta*gammalib.rad2deg))
    print("Distance model-pointing (lambda) .....: %.2f deg" % (lamda*gammalib.rad2deg))
    print("Omega_0 ..............................: %.2f deg" % (omega0*gammalib.rad2deg))
    print("PA of pointing from model ............: %.2f deg" % (pa_pnt*gammalib.rad2deg))
    print("Rho interval .........................: %.2f-%.2f deg" % (rho_min*gammalib.rad2deg,rho_max*gammalib.rad2deg))

    # Plot rhos
    for rho in rhos:

        # Print header
        print("Rho ..................................: %.2f deg" % (rho))

        # Compute some stuff
        rho_rad   = rho   * gammalib.deg2rad

        # Compute half length of arc that lies within Psf validity circle
        domega = 0.5 * cta_roi_arclength(rho_rad, zeta, delta_max)

        # Print domega
        print("  domega .............................: %.2f deg" % (domega*gammalib.rad2deg))

        # Continue only if arc length is positive
        if (domega > 0.0):

            # Precompute cosine and sine terms for azimuthal integration
            cos_rho = math.cos(rho_rad);
            sin_rho = math.sin(rho_rad);
            cos_psf = cos_rho * math.cos(zeta)
            sin_psf = sin_rho * math.sin(zeta)
            cos_ph  = cos_rho * math.cos(lamda)
            sin_ph  = sin_rho * math.sin(lamda)

            # Full containment within ellipse
            if (rho < semiminor):

                # Set Omega range
                omega_min = -domega
                omega_max = +domega

                # Print
                print("  Omega range ........................: %.2f - %.2f deg" % (omega_min*gammalib.rad2deg,omega_max*gammalib.rad2deg))

                # Plot segment
                n    = 100
                x    = []
                y    = []
                step = (omega_max - omega_min)/n
                for i in range(n+1):
                    radians = i*step + omega_min
                    cos_radians = math.sin(radians+obsOmega)
                    sin_radians = math.cos(radians+obsOmega)
                    x.append(rho * cos_radians)
                    y.append(rho * sin_radians)
                plt.plot(x, y, 'k-')

            # Partial containment within ellipse
            else:

                # Azimuth angle of intersection points between ellipse and
                # circle with radius of rho
                arg1        = 1.0 - (semiminor*semiminor) / (rho*rho)
                arg2        = 1.0 - (semiminor*semiminor) / (semimajor*semimajor)
                omega_width = math.acos(math.sqrt(arg1/arg2))

                # Print Omega width
                print("  omega_width ........................: %.2f deg" % (omega_width*gammalib.rad2deg))

                # 
                if (omega_width > 0.0):
                
                    # Compute Omega intervals. Recall that Omega is with respect
                    # to connecting line between model and observed photon. Make
                    # sure that omega_0 is within [-pi,pi] thus that the omega
                    # intervals are between [-2pi,2pi]
                    omega_0 = posang * gammalib.deg2rad - obsOmega
                    if (omega_0 > gammalib.pi):
                        omega_0 -= gammalib.pi
                    elif (omega_0 < -gammalib.pi):
                        omega_0 += gammalib.pi

                    # Compute intervals
                    omega1_min = omega_0    - omega_width
                    omega1_max = omega_0    + omega_width
                    omega2_min = omega1_min + gammalib.pi
                    omega2_max = omega1_max + gammalib.pi

                    print("  Omega range 1 (initial) ............: %.2f - %.2f deg" % (omega1_min*gammalib.rad2deg,omega1_max*gammalib.rad2deg))
                    print("  Omega range 2 (initial) ............: %.2f - %.2f deg" % (omega2_min*gammalib.rad2deg,omega2_max*gammalib.rad2deg))

                    # Limit each of the intervals to [omega_min, omega_max]
                    omega_min       = -domega
                    omega_max       = +domega
                    omega_min_plus  = omega_min + gammalib.twopi
                    omega_max_plus  = omega_max + gammalib.twopi
                    omega_min_minus = omega_min - gammalib.twopi
                    omega_max_minus = omega_max - gammalib.twopi

                    # If interval 1 overlaps then limit interval boundaries
                    if (omega1_max > omega_min and omega1_min < omega_max):
                        if (omega1_min < omega_min):
                            omega1_min = omega_min
                        if (omega1_max > omega_max):
                            omega1_max = omega_max
                    elif (omega1_max > omega_min_plus and omega1_min < omega_max_plus):
                        if (omega1_min < omega_min_plus):
                            omega1_min = omega_min_plus
                        if (omega1_max > omega_max_plus):
                            omega1_max = omega_max_plus
                    elif (omega1_max > omega_min_minus and omega1_min < omega_max_minus):
                        if (omega1_min < omega_min_minus):
                            omega1_min = omega_min_minus
                        if (omega1_max > omega_max_minus):
                            omega1_max = omega_max_minus
                    else:
                        omega1_min = 0.0 # Dummy for no-interval
                        omega1_max = 0.0

                    # If interval 2 overlaps then limit interval boundaries
                    if (omega2_max > omega_min and omega2_min < omega_max):
                        if (omega2_min < omega_min):
                            omega2_min = omega_min
                        if (omega2_max > omega_max):
                            omega2_max = omega_max
                    elif (omega2_max > omega_min_plus and omega2_min < omega_max_plus):
                        if (omega2_min < omega_min_plus):
                            omega2_min = omega_min_plus
                        if (omega2_max > omega_max_plus):
                            omega2_max = omega_max_plus
                    elif (omega2_max > omega_min_minus and omega2_min < omega_max_minus):
                        if (omega2_min < omega_min_minus):
                            omega2_min = omega_min_minus
                        if (omega2_max > omega_max_minus):
                            omega2_max = omega_max_minus
                    else:
                        omega2_min = 0.0 # Dummy for no-interval
                        omega2_max = 0.0

                    # Print results
                    print("  omega_0 ............................: %.2f deg" % (omega_0*gammalib.rad2deg))
                    print("  Omega range ........................: %.2f - %.2f deg" % (omega_min*gammalib.rad2deg,omega_max*gammalib.rad2deg))
                    print("  Omega range 1 ......................: %.2f - %.2f deg" % (omega1_min*gammalib.rad2deg,omega1_max*gammalib.rad2deg))
                    print("  Omega range 2 ......................: %.2f - %.2f deg" % (omega2_min*gammalib.rad2deg,omega2_max*gammalib.rad2deg))

                    # Draw Omega range lines and circle segments
                    if (omega1_min < omega1_max):
                        #n    = 10
                        #step = rho/n
                        #x    = []
                        #y    = []
                        #for i in range(n+1):
                        #    x.append(i*step * math.sin(omega1_min+obsOmega))
                        #    y.append(i*step * math.cos(omega1_min+obsOmega))
                        #plt.plot(x, y, 'k-.')
                        #x    = []
                        #y    = []
                        #for i in range(n+1):
                        #    x.append(i*step * math.sin(omega1_max+obsOmega))
                        #    y.append(i*step * math.cos(omega1_max+obsOmega))
                        #plt.plot(x, y, 'k-.')
                        n    = 100
                        x    = []
                        y    = []
                        step = (omega1_max - omega1_min)/n
                        for i in range(n+1):
                            radians = i*step + omega1_min
                            cos_radians = math.sin(radians+obsOmega)
                            sin_radians = math.cos(radians+obsOmega)
                            x.append(rho * cos_radians)
                            y.append(rho * sin_radians)
                        plt.plot(x, y, 'k-')

                    if (omega2_min < omega2_max):
                        #n    = 10
                        #step = rho/n
                        #x    = []
                        #y    = []
                        #for i in range(n+1):
                        #    x.append(i*step * math.sin(omega2_min+obsOmega))
                        #    y.append(i*step * math.cos(omega2_min+obsOmega))
                        #plt.plot(x, y, 'k-.')
                        #x    = []
                        #y    = []
                        #for i in range(n+1):
                        #    x.append(i*step * math.sin(omega2_max+obsOmega))
                        #    y.append(i*step * math.cos(omega2_max+obsOmega))
                        #plt.plot(x, y, 'k-.')
                        n    = 100
                        x    = []
                        y    = []
                        step = (omega2_max - omega2_min)/n
                        for i in range(n+1):
                            radians = i*step + omega2_min
                            cos_radians = math.sin(radians+obsOmega)
                            sin_radians = math.cos(radians+obsOmega)
                            x.append(rho * cos_radians)
                            y.append(rho * sin_radians)
                        plt.plot(x, y, 'k-')

        # Test line
        #x = []
        #y = []
        #for i in range(10):
        #    x.append(i*0.1 * math.sin(0.0+obsOmega))
        #    y.append(i*0.1 * math.cos(0.0+obsOmega))
        #plt.plot(x, y, 'k-')

        # Build circle
        #circle_x = []
        #circle_y = []
        #for angle in range(0, 360):
        #    radians     = angle * gammalib.deg2rad
        #    cos_radians = math.cos(radians)
        #    sin_radians = math.sin(radians)
        #    x = rho * cos_radians
        #    y = rho * sin_radians
        #    circle_x.append(x)
        #    circle_y.append(y)
        #plt.plot(circle_x, circle_y, 'k-.')

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
    print("****************************************")
    print("* Test elliptical response computation *")
    print("****************************************")

    # Check if matplotlib is available
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Matplotlib is not installed on your system. Abort.")
        sys.exit()

    # Set model parameters
    centre = gammalib.GSkyDir()
    centre.radec_deg(0.0, 0.0)
    semimajor = 1.0
    semiminor = 0.50
    posang    = -99.0

    # Set observed photon direction
    obsDir = gammalib.GSkyDir()
    obsDir.radec_deg(0.4, -0.3)

    # Set pointing direction
    pntDir = gammalib.GSkyDir()
    pntDir.radec_deg(0.1, -1.0)

    # Set rhos
    rhos = [i*0.02 for i in range(51)]

    # Show elliptical response geometry
    show_response(centre, semimajor, semiminor, posang, obsDir, pntDir, rhos)
    
