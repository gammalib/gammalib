#! /usr/bin/env python
# ==========================================================================
# Approximate radial disk function
# --------------------------------------------------------------------------
#
# Copyright (C) 2020 Juergen Knoedlseder
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
import math
import gammalib
import matplotlib.pyplot as plt


# ========================== #
# Compare radial model to MC #
# ========================== #
def show(theta0, k=1000.0, rad_max=1.0, nrad=1000):
    """
    Show smooth approximation of radial disk function.
    """
    # Get dummy energy and time
    energy = gammalib.GEnergy()
    time   = gammalib.GTime()
    
    # Generate theta vector
    theta_bin = rad_max/nrad
    theta     = [(i+0.5) * theta_bin for i in range(nrad)]

    # Generate model vector
    model = [1.0 / (1.0 + math.exp(-k * (theta0 - t))) for t in theta]
    disk  = [1.0 if t < theta0 else 0.0 for t in theta]

    # Create figure
    plt.figure(1)
    plt.title('Smooth approximation of radial disk function')

    # Plot data
    plt.plot(theta, model, 'r-')
    plt.plot(theta, disk, 'g-')

    # Set axes
    plt.xlabel('Offset (deg)')
    plt.ylabel('Arbitrary units')

    # Notify
    print('PLEASE CLOSE WINDOW TO CONTINUE ...')

    # Show plot
    plt.show()

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show radial model
    show(0.5)
