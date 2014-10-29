#! /usr/bin/env python
# ==========================================================================
# This script tests the radial response computation. The script requires
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
    

# ==================== #
# Show radial response #
# ==================== #
def show_response(centre, radius, obsDir, delta_max):
    """
    Show radial response
    """
    # Create figure
    plt.figure(1)
    plt.title("Radial response")

    # Compute some stuff
    rho_obs      = centre.dist_deg(obsDir)
    posangle_obs = centre.posang(obsDir)

    # Rho limits
    rho_min = 0.0
    rho_max = rho_obs + delta_max
    if (rho_obs > delta_max):
        rho_min = rho_obs - delta_max
    if rho_max > radius:
        rho_max = radius

    # Convert to radians
    delta_max *= gammalib.deg2rad
    rho_obs   *= gammalib.deg2rad

    # Setup rhos
    rhos = []
    for i in range(51):
        rho = i*0.02
        if (rho >= rho_min and rho <= rho_max):
            rhos.append(rho)

    # Plot radial model
    x = []
    y = []
    for angle in range(0, 361):
        radians     = angle * gammalib.deg2rad
        cos_radians = math.cos(radians)
        sin_radians = math.sin(radians)
        x.append(radius * cos_radians)
        y.append(radius * sin_radians)
    plt.plot(x, y, 'r-')
    plt.xlim([ 2.0*radius,-2.0*radius])
    plt.ylim([-2.0*radius, 2.0*radius])
    plt.plot([0.0], [0.0], 'ro')

    # Plot observed direction
    x_obs = obsDir.ra_deg()  - centre.ra_deg()
    y_obs = obsDir.dec_deg() - centre.dec_deg()
    plt.plot([x_obs], [y_obs], 'bo')

    # Draw connecting line model-photon'
    n    = 10
    step = rho_obs/n * gammalib.rad2deg
    x    = []
    y    = []
    for i in range(n+1):
        x.append(i*step * math.sin(posangle_obs))
        y.append(i*step * math.cos(posangle_obs))
    plt.plot(x, y, 'b-')

    # Plot Psf validity circle
    x = []
    y = []
    for angle in range(0, 361):
        radians     = angle * gammalib.deg2rad
        cos_radians = math.cos(radians)
        sin_radians = math.sin(radians)
        x.append(delta_max * cos_radians * gammalib.rad2deg + x_obs)
        y.append(delta_max * sin_radians * gammalib.rad2deg + y_obs)
    plt.plot(x, y, 'b-')

    # Plot rhos
    for rho in rhos:

        # Compute some stuff
        rho_rad = rho * gammalib.deg2rad

        # Compute half length of arc that lies within Psf validity circle
        domega = 0.5 * cta_roi_arclength(rho_rad, rho_obs, delta_max)

        # Continue only if arc length is positive
        if (domega > 0.0):

            # Set Omega range
            omega_min = -domega
            omega_max = +domega

            # Plot segment
            n    = 100
            x    = []
            y    = []
            step = (omega_max - omega_min)/n
            for i in range(n+1):
                radians = i*step + omega_min
                cos_radians = math.sin(radians+posangle_obs)
                sin_radians = math.cos(radians+posangle_obs)
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
    Test radial response computation
    """
    # Dump header
    print("")
    print("************************************")
    print("* Test radial response computation *")
    print("************************************")

    # Check if matplotlib is available
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Matplotlib is not installed on your system. Abort.")
        sys.exit()

    # Set Psf size
    delta_max = 0.6

    # Set model parameters
    centre = gammalib.GSkyDir()
    centre.radec_deg(0.0, 0.0)
    radius = 0.8

    # Set observed photon direction
    obsDir = gammalib.GSkyDir()
    obsDir.radec_deg(0.4, -0.3)
    #obsDir.radec_deg(-0.4, -0.3)
    #obsDir.radec_deg(0.4, 0.4)

    # Show elliptical response geometry
    show_response(centre, radius, obsDir, delta_max)
    
