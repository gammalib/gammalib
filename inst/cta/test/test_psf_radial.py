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
def show_response(centre, theta_max, obsDir, psf_max):
    """
    Show radial response
    """
    # Create figure
    plt.figure(1)
    plt.title("Radial response")

    # Compute some stuff
    delta_mod    = obsDir.dist_deg(centre)
    posangle_mod = obsDir.posang(centre)

    # Delta limits
    delta_min = 0.0
    delta_max = delta_mod + theta_max
    if (delta_mod > theta_max):
        delta_min = delta_mod - theta_max
    if delta_max > psf_max:
        delta_max = psf_max

    # Convert to radians
    delta_mod *= gammalib.deg2rad
    theta_max *= gammalib.deg2rad

    # Setup deltas
    deltas = []
    for i in range(51):
        delta = i*0.02
        if (delta >= delta_min and delta <= delta_max):
            deltas.append(delta)

    # Plot Psf validity circle
    x = []
    y = []
    for angle in range(0, 361):
        radians     = angle * gammalib.deg2rad
        cos_radians = math.cos(radians)
        sin_radians = math.sin(radians)
        x.append(psf_max * cos_radians)
        y.append(psf_max * sin_radians)
    plt.plot(x, y, 'b-')
    plt.xlim([ 2.0*psf_max,-2.0*psf_max])
    plt.ylim([-2.0*psf_max, 2.0*psf_max])
    plt.plot([0.0], [0.0], 'bo')

    # Plot model direction
    x_mod = centre.ra_deg()  - obsDir.ra_deg()
    y_mod = centre.dec_deg() - obsDir.dec_deg()
    plt.plot([x_mod], [y_mod], 'ro')

    # Plot radial model
    x = []
    y = []
    for angle in range(0, 361):
        radians     = angle * gammalib.deg2rad
        cos_radians = math.cos(radians)
        sin_radians = math.sin(radians)
        x.append(theta_max * cos_radians * gammalib.rad2deg + x_mod)
        y.append(theta_max * sin_radians * gammalib.rad2deg + y_mod)
    plt.plot(x, y, 'r-')

    # Plot deltas
    for delta in deltas:

        # Compute some stuff
        delta_rad = delta * gammalib.deg2rad

        # Compute half length of arc that lies within model
        dphi = 0.5 * cta_roi_arclength(delta_rad, delta_mod, theta_max)

        # Continue only if arc length is positive
        if (dphi > 0.0):

            # Set Phi range
            phi_min = -dphi
            phi_max = +dphi

            # Plot segment
            n    = 100
            x    = []
            y    = []
            step = (phi_max - phi_min)/n
            for i in range(n+1):
                radians = i*step + phi_min
                cos_radians = math.sin(radians+posangle_mod)
                sin_radians = math.cos(radians+posangle_mod)
                x.append(delta * cos_radians)
                y.append(delta * sin_radians)
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
    print("*******************************")
    print("* Test radial Psf computation *")
    print("*******************************")

    # Check if matplotlib is available
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Matplotlib is not installed on your system. Abort.")
        sys.exit()

    # Set Psf size
    psf_max = 1.0

    # Set observed photon direction
    obsDir = gammalib.GSkyDir()
    obsDir.radec_deg(0.0, 0.0)

    # Set model parameters
    centre = gammalib.GSkyDir()
    centre.radec_deg(0.0, 0.6)
    theta_max = 0.5

    # Show elliptical response geometry
    show_response(centre, theta_max, obsDir, psf_max)
    
