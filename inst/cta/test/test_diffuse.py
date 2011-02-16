#! /usr/bin/env python

from gammalib import *


# ====================== #
# Set point source model #
# ====================== #
def ptsrc_model(ra=0.0, dec=0.0):
	"""
	Set shell model.
	"""
	# Set shell centre
	pos = GSkyDir()
	pos.radec_deg(ra, dec)
	
	# Set spatial model
	spatial = GModelSpatialPtsrc(pos)
	
	# Set spectral model
	spectral = GModelSpectralPlaw(1.0, -2.0)
	
	# Set sky model
	model = GModelPointSource(spatial, spectral)
	
	# Optionally show model
	#print model
	
	# Return model
	return model


# =============== #
# Set shell model #
# =============== #
def shell_model(ra=0.0, dec=0.0, radius=1.0, width=0.3):
	"""
	Set shell model.
	"""
	# Set shell centre
	center = GSkyDir()
	center.radec_deg(ra, dec)
	
	# Set spatial model
	spatial = GModelSpatialShell(center, radius, width, False)
	
	# Set spectral model
	spectral = GModelSpectralPlaw(1.0, -2.0)
	
	# Set sky model
	model = GModelDiffuseSource(spatial, spectral)
	
	# Optionally show model
	#print model
	
	# Return model
	return model


# ========================== #
# Set binned CTA observation #
# ========================== #
def observation(ra=0.0, dec=0.0, binsz=0.02, npix=100, ebins=10):
	"""
	Set binned CTA observation.
	"""
	# Allocate observation
	obs = GCTAObservation()
	
	# Set response
	obs.response("kb_E_50h_v3", "../caldb")
	
	# Set pointing
	dir = GSkyDir()
	pnt = GCTAPointing()
	dir.radec_deg(ra, dec)
	pnt.dir(dir);
	obs.pointing(pnt)
	
	# Set
	ebounds = GEbounds()
	emin    = GEnergy()
	emax    = GEnergy()
	emin.TeV(0.1)
	emax.TeV(100.0)
	ebounds.setlog(emin, emax, ebins)
	gti     = GGti()
	tmin    = GTime()
	tmax    = GTime()
	tmin.met(0.0)
	tmax.met(1800.0)
	gti.append(tmin, tmax)
	map     = GSkymap("CAR", "CEL", ra, dec, -binsz, binsz, npix, npix, ebins)
	cube    = GCTAEventCube(map, ebounds, gti)
	obs.events(cube)
	
	# Optionally show observation
	#print obs
	
	# Return observation
	return obs


# ==================== #
# Test IRF computation #
# ==================== #
def test_irf(model, filename="cntmap.fits"):
	"""
	Test IRF.
	"""
	# Set CTA observation
	obs = observation()
	src = GSkyDir()
	src.radec_deg(0.0,0.0)
	
	# Loop over all bins
	for bin in obs.events():
		
		# Cast to CTA bin
		bin = cast_GCTAEventBin(bin)
		
		# Set bin energy and time as source energy and time (no dispersion)
		srcEng  = bin.energy()
		srcTime = bin.time()
		
		# Compute IRF
		irf = obs.response().irf(bin, model, srcEng, srcTime, obs) #*bin.size()
		
		# Set bin
		bin.counts(irf)
		
		# Optionally print result
		#print bin.energy(), bin.dir(), bin.time(), irf
		
		# Compute distance
		#if irf > 0.0:
		#	distance = bin.dir().dist(src)
		#	sigma    = obs.response().psf_dummy_sigma(bin.energy().log10TeV())
		#	print distance/sigma, irf
	
	# Save observation
	obs.save(filename, True)
	
	# Return
	return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
	"""
	Test diffuse models.
	"""
	# Dump header
	print
	print "***********************"
	print "* Test diffuse models *"
	print "***********************"
	
	# Set shell model
	#model = ptsrc_model()
	model = shell_model()
	
	# Test IRF
	test_irf(model)

