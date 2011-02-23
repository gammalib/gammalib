#! /usr/bin/env python
from gammalib import *

def test_radial_model(name, model, clobber=True):
	"""Make an image of a spatial model 
	and check that it integrates to 1"""
	# Make a 2D FITS image
	ra, dec, binsz, npix = 0, 0, 0.01, 200
	image = GSkymap("CAR", "CEL", ra, dec, -binsz, binsz, npix, npix, 1)	

	# Fill the image	
	for pix in range(image.npix()):
		dir = image.pix2dir(pix)
		theta = center.dist(dir)
		image[pix, 0] = model.eval(theta);
	
	# Write it to file
	filename = name + '.fits'
	print 'Writing', filename
	image.save(filename, clobber)

	# Check that it integrates to 1
	integral = 0
	for pix in range(image.npix()):
		integral += image.omega(pix) * image[pix]
	print 'integral = %g (should be 1)' % integral 
	print 'integral error = %g (should be 0)' % (integral - 1)
	print

if __name__ == '__main__':
	
	center = GSkyDir()
	center.radec_deg(0, 0)

	test_radial_model(name='gauss', model=GModelRadialGauss(center, 0.3))
	test_radial_model(name='disk', model=GModelRadialDisk(center, 0.8))
	test_radial_model(name='shell', model=GModelRadialShell(center, 0.5, 0.1))
