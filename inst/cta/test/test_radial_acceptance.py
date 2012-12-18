#! /usr/bin/env python
# ==========================================================================
# This script displays the radial acceptance model that may be used to model
# the CTA radial acceptance.
#
# Requires:
# - matplotlib
#
# Copyright (C) 2011-2012 Juergen Knoedlseder
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
import matplotlib.pyplot as plt
from gammalib import *
from math import *


# ========== #
# Show model #
# ========== #
def show_model(model):
    """
    Show radial acceptance model using matplotlib.
    """
    # Create angular axis (from 0 to 4 deg)
    thetas = [i * 0.05 for i in range(80)]

    # Get model values
    values      = [model.eval(theta) for theta in thetas]
    values_grad = [model.eval_gradients(theta) for theta in thetas]

    # Create figure
    plt.figure(1)
    plt.title("Radial acceptance model (" + model.type() + ")")

    # Plot data
    plt.plot(thetas, values, 'r-')
    plt.plot(thetas, values_grad, 'ro')

    # Set axes
    plt.xlabel("Offset angle (deg)")
    plt.ylabel("Function value")

    # Show plot
    plt.show()

    # Return
    return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    Show radial acceptance models.
    """
    # Dump header
    print
    print "*************************************"
    print "* Show CTA radial acceptance models *"
    print "*************************************"

    # Display various models
    show_model(GCTAModelRadialGauss(3.0))
    show_model(GCTAModelRadialProfile(1.5, 3.0, 5.0))
    show_model(GCTAModelRadialPolynom([1.0, -0.1239176, +0.9751791,
                                       -3.0584577, +2.9089535, -1.3535372,
                                       +0.3413752, -0.0449642, +0.0024321]))
