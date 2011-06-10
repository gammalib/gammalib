#! /usr/bin/env python
# ===========================================================================================#
# This script makes a model map from an XML file.
# ===========================================================================================#
from gammalib import *


# ============================ #
# Radial model testing routine #
# ============================ #
def make_model_for_cntmap(cntmap, modname, xmlname, irf, caldb, clobber=True):
	"""
	Make a model counts map from an XML model.
	"""
	# Allocate empty CTA observation
	obs = GCTAObservation()
	
	# Load counts map into CTA observation
	obs.load_binned(cntmap)
	
	# Set response
	obs.response(irf, caldb)

	# Load models from XML file
	models = GModels(xmlname)
	
	# Loop over all bins in counts map
	events = cast_GCTAEventCube(obs.events())
	for event in events:
		model = models.eval(event, obs) * event.size()
		event.counts(model)
	
	# Save CTA observation
	obs.save(modname, clobber)

	# Return
	return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
	"""
	Example illustrating how to create a model map for a counts map.
	"""
	# Dump header
	print
	print "***************************************************"
	print "* Create model map for counts map using XML model *"
	print "***************************************************"
	print "... please wait for a few seconds"

	# Set parameters
	irf     = "kb_E_50h_v3"
	caldb   = "../caldb"
	xmlname = "data/crab.xml"
	cntmap  = "data/crab_cntmap.fits"
	modname = "model.fits"
	
    # Make the model
	make_model_for_cntmap(cntmap, modname, xmlname, irf, caldb)
	
