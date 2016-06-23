#! /usr/bin/env python
# ==========================================================================
# This script shows the shell model radial distribution.
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
# Check if matplotlib is available
try:
    import matplotlib.pyplot as plt
    has_matplotlib = True
except ImportError:
    print("Matplotlib is not (correctly) installed on your system.")
    has_matplotlib = False


# ================ #
# Compute function #
# ================ #
def function(theta, radius=0.30, width=0.10):
    """
    """
    # Initialise
    value          = 0.0
    sin_theta      = math.sin(theta*gammalib.deg2rad)
    sin_theta2     = sin_theta*sin_theta
    sin_theta_in   = math.sin(radius*gammalib.deg2rad)
    sin_theta2_in  = sin_theta_in*sin_theta_in
    sin_theta_out  = math.sin((radius+width)*gammalib.deg2rad)
    sin_theta2_out = sin_theta_out*sin_theta_out

    # Compute function value
    if (sin_theta2 < sin_theta2_out):
        value = math.sqrt(sin_theta2_out - sin_theta2)
        if (sin_theta2 < sin_theta2_in):
            value -= math.sqrt(sin_theta2_in - sin_theta2)

    # Return value
    return value


# ============= #
# Show function #
# ============= #
def show(theta_max=0.5, theta_bins=100):
    """
    """

    # Continue only if matplotlib is available
    if has_matplotlib:

        # Create figure
        plt.figure(1)
        plt.title("Shell function")

        # Create theta axis
        theta  = []
        dtheta = theta_max/theta_bins
        for i in range(theta_bins):
            theta.append((i)*dtheta)

        # Get function
        value = []
        for i in range(theta_bins):
            value.append(function(theta[i]))

        # Plot function
        plt.plot(theta, value, 'r-')
        
        # Set axes
        plt.xlabel("Theta angle (degrees)")
        plt.ylabel("Function")

        # Create figure
        plt.figure(2)
        plt.title("MC function")

        # Get functions
        value    = []
        envelope = []
        for i in range(theta_bins):
            value.append(function(theta[i])*math.sin(theta[i]*gammalib.deg2rad))
            envelope.append(math.sin(0.4*gammalib.deg2rad)*math.sin(0.4*gammalib.deg2rad))

        # Plot function
        plt.plot(theta, value, 'r-')
        plt.plot(theta, envelope, 'b-')
        
        # Set axes
        plt.xlabel("Theta angle (degrees)")
        plt.ylabel("Function*sin(theta)")

        # Notify
        print("PLEASE CLOSE WINDOW TO CONTINUE ...")

        # Show plot
        plt.show()

    # Return
    return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    Show Shell function profile
    """
    # Dump header
    print("")
    print("*******************************")
    print("* Show Shell function profile *")
    print("*******************************")

    # Show
    show()
