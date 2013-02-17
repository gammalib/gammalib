#! /usr/bin/env python
# ==========================================================================
# This script makes an image from a radial model and verifies that the model
# is correctly normalized.
# ==========================================================================
from gammalib import *


# ============================ #
# Radial model testing routine #
# ============================ #
def test_radial_model(name, model, clobber=True):
    """Make an image of a spatial model
    and check that it integrates to 1"""
    # Make a 2D FITS image (cartesian projection, celestial coordinates)
    ra, dec, binsz, npix = 0, 0, 0.01, 300
    image = GSkymap("CAR", "CEL", ra, dec, -binsz, binsz, npix, npix, 1)

    # Fill the image
    for pix in range(image.npix()):
        dir = image.pix2dir(pix)
        theta = center.dist(dir)
        image[pix] = model.eval(theta)

    # Write it to file
    filename = name + '.fits'
    print 'Writing', filename
    image.save(filename, clobber)

    # Check that it integrates to 1
    integral = 0.0
    for pix in range(image.npix()):
        integral += image.omega(pix) * image[pix]
    print 'integral = %g (should be 1)' % integral
    print 'integral error = %g (should be 0)' % (integral - 1)
    print

# ================ #
# Main entry point #
# ================ #
if __name__ == '__main__':
    """
    Test radial model integration.
    """
    # Set model location / centre
    center = GSkyDir()
    center.radec_deg(0.3, 0.1)

    # Test models
    test_radial_model(name='gauss', model=GModelSpatialRadialGauss(center, 0.3))
    test_radial_model(name='disk',  model=GModelSpatialRadialDisk(center, 0.8))
    test_radial_model(name='shell', model=GModelSpatialRadialShell(center, 0.5, 0.1))
