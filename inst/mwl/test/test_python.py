#! /usr/bin/env python

from gammalib import *
from math import *
import os


# ================================== #
# Extract data points from FITS file #
# ================================== #
def extract_data(filename):
	"""
	Extract data from FITS file.
	"""
	# Set plot styles
	ecolor = ['r',  'b',  'g']
	styles = ['ro', 'bo', 'go']
	
	# Open FITS file
	fits = GFits(filename)
	
	# Loop over all table HDUs
	k = 0
	for extno in range(fits.size()):
		if fits.hdu(extno).exttype() == GFitsHDU.HT_ASCII_TABLE or \
		   fits.hdu(extno).exttype() == GFitsHDU.HT_BIN_TABLE:
			
			# Load spectrum
			spectrum = GMWLSpectrum()
			spectrum.load(filename, extno)
			
			# Add spectral points
			x  = []
			y  = []
			dy = []
			for i in range(spectrum.size()):
				xval = spectrum[i].energy().MeV()
				yval = spectrum[i].flux()
				yerr = spectrum[i].flux_err()
				conv = xval*xval*1.6021765e-6
				x.append(xval)
				y.append(yval*conv)
				dy.append(yerr*conv)
			plt.loglog(x, y, styles[k])
			plt.errorbar(x, y, dy, fmt=None, ecolor=ecolor[k])

			# Increment style counter
			k = k + 1
	
	# Add axes
	plt.xlabel("Energy (MeV)")
	plt.ylabel("Flux (erg/cm2/s)")
	
	# Close FITS file
	fits.close()


# ============ #
# Fit spectrum #
# ============ #
def fit_spectrum(filename, xmlname):
	"""
	Fit spectrum.
	"""
	# Load observations
	obs     = GObservations()
	comptel = GMWLObservation(filename, "COMPTEL");
	lat     = GMWLObservation(filename, "LAT");
	hess    = GMWLObservation(filename, "HESS");
	obs.append(comptel)
	obs.append(lat)
	obs.append(hess)
	
	# Setup model
	models = GModels(xmlname)
	obs.models(models)
	
	# Perform optimization
	opt = GOptimizerLM()
	opt.max_iter(1000)
	obs.optimize(opt)
	print obs
	
	# Get model values
	x = np.power(10., np.arange(-1., 8., 0.1))
	y = []
	d = GSkyDir()
	d.radec_deg(83.6331, 22.0145)
	e = GEnergy()
	t = GTime()
	for energy in x:
		e.MeV(energy)
		yval = 0.0
		for model in obs.models():
			if model.type() == "PointSource":
			#crab = cast_GModelSky(model)
			#if crab != None:
				yval += model.value(d, e, t) * energy*energy*1.6021765e-6
		y.append(yval)
	plt.plot(x, y)
    

#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
	"""
	Perform testing.
	"""
	# Dump result
	print
	print "****************************"
	print "* Python interface testing *"
	print "****************************"
	
	# Only proceed if matplotlib is available
	try:
		# Import matplotlib
		import numpy as np
		import matplotlib.pyplot as plt

		# Open plot
		plt.figure(1)
		plt.title('GammaLib fit of the Crab PWN')

		# Extract data from FITS file
		extract_data("data/crab_mwl.fits")

		# Fit spectrum
		fit_spectrum("data/crab_mwl.fits", "data/crab_mwl2.xml")

		# Show spectrum
		plt.show()
		
	except ImportError:
		print "Matplotlib is not (correctly) installed on your system."
