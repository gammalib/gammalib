#! /usr/bin/env python
# ==========================================================================
# This script tests the Gaussian gradient.
#
# Requires:
# - matplotlib (optional)
#
# Copyright (C) 2012 Juergen Knoedlseder
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
from gammalib import *
from math import *


# ====================== #
# Test Gaussian gradient #
# ====================== #
def gauss(dir, sigma):
    """
    Test Gaussian gradient.
    """
    # Compute derivative
    dh = 0.0001
    s  = GSkyDir()
    e  = GEnergy()
    t  = GTime()
    f1 = GModelSpatialRadialGauss(s, sigma)
    f2 = GModelSpatialRadialGauss(s, sigma + dh)
    v1 = f1.eval_gradients(GPhoton(dir, e, t))
    g1 = f1[2].gradient()
    v2 = f2.eval_gradients(GPhoton(dir, e, t))
    g2 = f2[2].gradient()
    g  = (v2 - v1) / dh

    # Print result
    #print v1, v2
    #print g1, g2, g

    # Return
    return g


# ============= #
# Show Gaussian #
# ============= #
def show_gaussian(sigma):
    """
    Show Gaussian using matplotlib (if available).
    """
    # Only proceed if matplotlib is available
    try:
        # Import matplotlib
        import matplotlib.pyplot as plt

        # Create figure
        plt.figure(1)
        plt.title("Gaussian (sigma=" + str(sigma) + " deg)")

        # Setup gaussian
        skydir = GSkyDir()
        gauss  = GModelSpatialRadialGauss(skydir, sigma)

        # Create angular axis
        theta = [i * sigma * 0.05 for i in range(50)]

        # Extract function
        f_gauss    = []
        f_expected = []
        sigma_rad  = sigma * (pi / 180.0)
        norm       = 1.0 / (2.0 * pi * sigma_rad * sigma_rad)
        eng        = GEnergy()
        time       = GTime()
        for t in theta:
            s = GSkyDir()
            s.radec_deg(0.0, t)
            f = gauss.eval(GPhoton(s, eng, time))
            e = norm * exp(-0.5 * t * t / sigma / sigma)
            f_gauss.append(f)
            f_expected.append(e)

        # Plot data
        plt.plot(theta, f_gauss, 'r-')
        plt.plot(theta, f_expected, 'b.')

        # Set axes
        plt.xlabel("Separation (deg)")
        plt.ylabel("Function value")

        # Show plot
        plt.show()

    except ImportError:
        print("Matplotlib is not (correctly) installed on your system.")

    # Return
    return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    Test Gaussian gradient.
    """
    # Dump header
    print("")
    print("**************************")
    print("* Test Gaussian gradient *")
    print("**************************")

    # Show Gaussian
    show_gaussian(3.0)
