#! /usr/bin/env python
# ==========================================================================
# This script displays the CTA PSF using a response table.
#
# Requires:
# - matplotlib
# - numpy
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
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


# ====================== #
# Display response table #
# ====================== #
def display_response_table(index=1):
    """
    Display CTA response table.
    
    Keywords:
     index - Response parameter to display (starting from 0)
    """
    # Load PSF into response table
    fits  = GFits('../caldb/data/cta/e/bcf/000001/irf_test.fits')
    table = fits.table('POINT SPREAD FUNCTION')
    psf   = GCTAResponseTable(table)
    print(psf)

    # Set limits for 2D plot
    emin     = 0.05
    emax     = 80.0
    thetamin = 0.0
    thetamax = 2.5
    neng     = 100
    ntheta   = 100
    logemin  = log10(emin)
    logemax  = log10(emax)
    dlogeng  = (logemax - logemin) / (neng - 1)
    dtheta   = (thetamax - thetamin) / (ntheta - 1)

    # Fill image. Note that X and Y are inversed with respect to
    # traditional orientation
    image = np.zeros((ntheta, neng))
    for ieng in range(neng):
        eng = pow(10.0, logemin + ieng * dlogeng)
        for itheta in range(ntheta):
            theta               = thetamin + itheta * dtheta
            pars                = psf(eng, theta)
            image[itheta, ieng] = pars[index]

    def format_coord(x, y):
        col = int(x + 0.5)
        row = int(y + 0.5)
        if col >= 0 and col < ntheta and row >= 0 and row < neng:
            z = image[row, col]
            return 'x=%1.4f, y=%1.4f, z=%1.4f' % (x, y, z)
        else:
            return 'x=%1.4f, y=%1.4f' % (x, y)

    # Show image
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.imshow(image, cmap=cm.jet, interpolation='nearest', origin='lower')
    ax.format_coord = format_coord
    plt.show()

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    Perform testing.
    """
    # Dump result
    print("")
    print("******************************")
    print("* Display CTA response table *")
    print("******************************")

    # Display response table
    display_response_table()
