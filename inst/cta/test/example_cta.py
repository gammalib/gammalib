#! /usr/bin/env python

from gammalib import *
from math import *


# =================== #
# CTA binned analysis #
# =================== #
def binned_analysis(xmlname, cntmap):
	"""
	CTA binned analysis.
	"""
	# Load CTA observation
	obs = GCTAObservation()
	obs.load_binned(cntmap)

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
	print "*********************"
	print "* Test CTA analysis *"
	print "*********************"

    # Test
	binned_analysis("data/crab.xml", "data/crab_cntmap.fits")
