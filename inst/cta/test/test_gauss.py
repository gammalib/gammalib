#! /usr/bin/env python

from gammalib import *
from math import *


# ====================== #
# Test Gaussian gradient #
# ====================== #
def gauss(dir, sigma):
	"""
	Test Gaussian gradient.
	"""
	# Compute derivative
	dh = 0.0001
	s  = GSkyDir()
	f1 = GModelSpatialGauss(s, sigma)
	f2 = GModelSpatialGauss(s, sigma+dh)
	v1 = f1.eval_gradients(dir)
	g1 = f1[2].gradient()
	v2 = f2.eval_gradients(dir)
	g2 = f2[2].gradient()
	g  = (v2-v1)/dh

	# Print result
	print v1, v2
	print g1, g2, g

	# Return
	return g


# ============= #
# Show Gaussian #
# ============= #
def show_gaussian(sigma):
	"""
	Show Gaussian using matplotlib (if available).
	"""
	# Only proceed if matplotlib is available
	try:
		# Import matplotlib
		import matplotlib.pyplot as plt
		
		# Create figure
		plt.figure(1)
		plt.title("Gaussian (sigma="+str(sigma)+" deg)")
		
		# Setup gaussian
		skydir = GSkyDir()
		gauss  = GModelSpatialGauss(skydir, sigma)
		
		# Create angular axis
		theta = [i*sigma*0.05 for i in range(50)]
		
		# Extract function
		f_gauss    = []
		f_expected = []
		for t in theta:
			s = GSkyDir()
			s.radec_deg(0.0, t)
			f = gauss.eval_gradients(s)
			e = exp(-0.5*t*t/sigma/sigma)
			f_gauss.append(f)
			f_expected.append(e)
		
		# Plot data
		plt.plot(theta, f_gauss, 'r-')
		plt.plot(theta, f_expected, 'b.')
		
		# Set axes
		plt.xlabel("Separation (deg)")
		plt.ylabel("Function value")
		
		# Show plot
		plt.show()

	except:
		print "Matplotlib is not (correctly) installed on your system. No data are shown."

	# Return
	return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
	"""
	Test Gaussian gradient.
	"""
	# Dump header
	print
	print "**************************"
	print "* Test Gaussian gradient *"
	print "**************************"

    # Test
	#dir = GSkyDir()
	#dir.radec_deg(0.0, 1.0)
	#gauss(dir, 10.0)
	#
	show_gaussian(3.0)
