#! /usr/bin/env python

from gammalib import *
from math import *
import os
import sys


# ======================= #
# Test CTA response table #
# ======================= #
def cta_response_table(index=1):
    """
    """
    # Import modules for plotting
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    #
    t = GCTAResponseTable()
    print t

    #
    fits = GFits('/Users/jurgen/cta/dc1/psf/psf.fits')
    table = fits.table('POINT SPREAD FUNCTION')
    psf = GCTAResponseTable(table)
    print psf

    # Set limits for 2D plot
    emin = 0.05
    emax = 80.0
    thetamin = 0.0
    thetamax = 2.5
    neng = 100
    ntheta = 100
    logemin = log10(emin)
    logemax = log10(emax)
    dlogeng = (logemax - logemin) / (neng - 1)
    dtheta = (thetamax - thetamin) / (ntheta - 1)

    # Fill image. Note that X and Y are inversed with respect to
    # traditional orientation
    image = np.zeros((ntheta, neng))
    for ieng in range(neng):
        eng = pow(10.0, logemin + ieng * dlogeng)
        for itheta in range(ntheta):
            theta = thetamin + itheta * dtheta
            pars = psf(eng, theta)
            # print eng, theta, pars[index]
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
    ax = fig.add_subplot(111)
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
    print
    print "********************************"
    print "* CTA Python interface testing *"
    print "********************************"

    # Initialise success counter
    tests = 0
    success = 0

    # Perform tests
    cta_response_table()
