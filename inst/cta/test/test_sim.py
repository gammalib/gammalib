#! /usr/bin/env python

from gammalib import *
from math import *
import os
import numpy as np
import matplotlib.pyplot as plt


# ======================== #
# Simulate CTA observation #
# ======================== #
def simulate(xmlname):
	"""
	Simulate CTA observation.
	"""
    # Define MC parameters
	area   = 1.0
	dir    = GSkyDir()
	radius = 10.0
	emin   = GEnergy()
	emax   = GEnergy()
	tmin   = GTime()
	tmax   = GTime()
	emin.TeV(0.1)
	emax.TeV(100.0)
	tmin.met(0.0)
	tmax.met(1000.0)
	
	# Allocate random number generator
	ran = GRan()
	
	# Load models
	models = GModels(xmlname)
	print models
	
	# Simulate photons
	photons = models.mc(area, dir, radius, emin, emax, tmin, tmax, ran)
	
	# Print photons
	print str(len(photons))+" photon simulated."
	#for photon in photons:
	#	print photon

	# Return
	return

#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
	"""
	Simulate CTA observation.
	"""
	# Dump header
	print
	print "****************************"
	print "* Simulate CTA observation *"
	print "****************************"
	
    # Perform simulation
	simulate("data/source1.xml")
