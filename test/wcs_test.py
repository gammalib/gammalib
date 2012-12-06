#! /usr/bin/env python
# ===========================================================================================#
# This script creates WCS data structures using pywcs. It only works if pywcs, numpy and
# pyfits are installed. The main purpose of this script is to test the gammalib wcs interface
# against the wcslib interface.
# ===========================================================================================#
import numpy
import pywcs
import pyfits
import sys


# =================== #
# Test CAR projection #
# =================== #
def test_car(ra=83.63, dec=22.01):
    """
    Test CAR projection.
    """
    # Create new WCS object
    wcs = pywcs.WCS(naxis=2)

    # Set up projection
    wcs.wcs.crpix = [100.5, 100.5]
    wcs.wcs.cdelt = numpy.array([0.02, 0.02])
    wcs.wcs.crval = [ra, dec]
    wcs.wcs.ctype = ["RA---CAR", "DEC--CAR"]
    # wcs.wcs.ctype = ["GLON-CAR", "GLAT-CAR"]
    print wcs.to_header()

    # Debug
    print wcs.wcs.name
    # wcs.wcs.print_contents()

    # Write out the WCS object as FITS header
    header = wcs.to_header()
    # print header
    if dec >= 0.0:
        print "LONPOLE = 0"
        print "LATPOLE =", 90.0 - dec
    else:
        print "LONPOLE = 180"
        print "LATPOLE =", 90.0 + dec

    # Some pixel coordinates
    # pixcrd = numpy.array([[99.5,99.5],[100.5,100.5]], numpy.float_)
    pixcrd = numpy.array([[100.5, 100.5]], numpy.float_)
    sky = wcs.wcs_pix2sky(pixcrd, 0)
    print sky
    res = wcs.wcs.p2s(pixcrd, 0)
    print res
    res = wcs.wcs.s2p(numpy.array([[83.7163, 22.09]], numpy.float_), 0)
    print res


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    Perform testing.
    """
    # Dump result
    print
    print "**********************"
    print "* Python WCS testing *"
    print "**********************"

    # Initialise success counter
    tests = 0
    success = 0

    # Perform tests
    test_car()
