#! /usr/bin/env python
from gammalib import *

def test_GModelSpatialShell(ra=0, dec=0, radius=0.8, width=0.2,
							npix=300, small_angle=False,
							filename='shell.fits', clobber=True):
	"""Make a 2D image of the shell model 
	and check that it integrates to 1"""
	# Compute binsize
	binsz = 2.1*radius / npix

	# Set up model
	center = GSkyDir()
	center.radec_deg(ra, dec)

	model = GModelSpatialShell(center, radius, width, small_angle)

	# Make a 2D FITS image
	image = GSkymap("CAR", "CEL", ra, dec, -binsz, binsz, npix, npix, 1)
	for pix in range(image.npix()):
		dir = image.pix2dir(pix)
		image[pix, 0] = model.eval(dir);
	print 'Writing', filename
	image.save(filename, clobber)

	# Check that it integrates to 1
	integral = 0
	for pix in range(image.npix()):
		integral += image.omega(pix) * image[pix]
	print 'integral = %g (should be 1)' % integral 
	print 'integral error = %g (should be 0)' % (integral - 1)

if __name__ == '__main__':
	test_GModelSpatialShell(filename='small_angle_shell.fits', radius=20, small_angle=True)
	test_GModelSpatialShell(filename='large_angle_shell.fits', radius=20, small_angle=False)
