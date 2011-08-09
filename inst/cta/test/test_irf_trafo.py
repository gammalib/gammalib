#! /usr/bin/env python
# ===========================================================================================#
# Test coordinate transformations for IRF computation.
# ===========================================================================================#
from gammalib import *
from math import *


# ====================== #
# Compute Omega0 #
# ====================== #
def omega0():
	"""
	...
	"""
	# Set pointing direction
	pointing = GSkyDir()
	pointing.radec_deg(10.0, 10.0)

	# Set model centre
	model = GSkyDir()
	model.radec_deg(12.0, 10.0)

	# Set measured photon direction
	p = GSkyDir()
	p.radec_deg(13.0, 9.0)
	
	# Determine angular distance between measured photon direction and model centre
	zeta = model.dist(p)
	
	# Determine angular distance between measured photon direction and pointing direction
	eta = pointing.dist(p)
	
	# Determine angular distance between model centre and pointing direction
	lamda = model.dist(pointing)
	
	# Compute Omega0
	omega0 = 0.0
	omega2 = 0.0
	denom  = sin(lamda) * sin(zeta)
	if denom != 0.0:
		arg    = (cos(eta) - cos(lamda) * cos(zeta))/denom
		arg2   = cos(eta)/denom - 1./(tan(lamda) * tan(zeta))
		omega0 = acos(arg)
		omega2 = acos(arg2)
	
	# Print results
	print "Omega0 ...: ", omega0*180.0/pi
	print "Omega0 ...: ", omega2*180.0/pi
	print "Zeta .....: ", zeta
	print "Eta ......: ", eta
	print "Lambda ...: ", lamda
	
	# Return
	return



#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
	"""
	...
	"""
	# ...
	omega0()
