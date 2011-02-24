#! /usr/bin/env python

from gammalib import *
from math import *


# =================== #
# CTA binned analysis #
# =================== #
def binned_analysis(model, cntmap, irf, caldb):
	"""
	CTA binned analysis.
	"""
	# Allocate empty observation container
	obs = GObservations()
	
	# Allocate empty CTA observation
	cta_obs = GCTAObservation()
	
	# Load counts map into CTA observation
	cta_obs.load_binned(cntmap)
	
	# Specify response for CTA observation
	cta_obs.response(irf, caldb)
	
	# Append CTA observation to observation container
	obs.append(cta_obs)

	# Load model to describe the data from XML file
	obs.models(model)
	
	# Allocate Levenberg-Marquardt optimizer
	opt = GOptimizerLM()

	# Optimize model parameters
	obs.optimize(opt)

	# Get maximum likelihood value (just for fun)
	logL = -(opt.value())

	# Get a copy of the model fitting results. We want a copy
	# here as the models are part of the observation container
	# "obs", and the container goes out of scope once the function
	# is left (and thus the models would also get out of scope)
	models = obs.models().copy()
	
	# Print optimizer results
	print opt
	
	# Return models
	return models


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
	"""
	Example illustrating a binned data analysis.
	"""
	# Dump header
	print
	print "*********************************************************"
	print "* Perform binned maximum likelihood fitting of CTA data *"
	print "*********************************************************"
	print "... please wait for a few seconds"

	# Set parameters
	irf    = "kb_E_50h_v3"
	caldb  = "../caldb"
	model  = "data/crab.xml"
	cntmap = "data/crab_cntmap.fits"
	
    # Perform binned analysis
	result = binned_analysis(model, cntmap, irf, caldb)
	
	# Print model results
	print result
	
