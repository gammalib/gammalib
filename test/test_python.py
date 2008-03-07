#! /usr/bin/env python

from gammalib import *


#================#
# Test FITS file #
#================#
def test_fits():
	"""
	Test GammaLib GFits interface.
	"""
	# Allocate empty FITS file
	fits = GFits()
	
	# Print empty FITS file
	#print fits

	# Open FITS file
	fits.open("test_python.fits")
	
	# Print real FITS file
	#print fits
	
	# Print number of HDUs
	print "Number of HDUs: " + str(fits.num_hdus())
	
	# Allocate HDU and append it to FITS file
	hdu = GFitsHDU()
	fits.append_hdu(hdu)
	#print hdu
	#print fits
	
	# Save FITS file
	fits.save()
	
	# Save into another FITS file
	fits.saveto("test_python_2.fits")


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
	"""
	Perform testing.
	"""
	test_fits()
