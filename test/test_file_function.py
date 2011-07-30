#! /usr/bin/env python
# ==========================================================================
# This script tests the GModelSpectralFunc file function spectral model.
#
# Copyright (C) 2011 Jurgen Knodlseder
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ==========================================================================
from gammalib import *
import math


# ================================= #
# Generate file function from model #
# ================================= #
def make_file_function(model, emin=0.1e6, emax=100.0e6, nodes=10):
	"""
	Generate file function from model.
	"""
	# Set filename
	filename = "file_function.dat"
	
	# Generate file function CSV file
	dloge = (math.log10(emax)-math.log10(emin))/(nodes-1)
	file = open(filename, "w")
	for i in range(nodes):
		energy = math.pow(10.0, math.log10(emin) + i*dloge)
		eng    = GEnergy()
		eng.MeV(energy)
		value  = model.eval(eng)
		file.write(str(energy)+" "+str(value)+"\n")
	file.close()
	
	# Generate file function
	file_function = GModelSpectralFunc(filename)
	
	# Return
	return file_function


# ========================== #
# Compare photon flux values #
# ========================== #
def compare_photon_flux(model, file_function, emin=8.7, emax=1123.0, steps=13):
	"""
	Compare photon fluxes
	"""
	# Set accuracy
	eps       = 1.0e-8
	n_violate = 0
	
	# Allocate boundaries
	e_min = GEnergy()
	e_max = GEnergy()
	
	# Ramp up emax from below first node
	e_min.MeV(emin)
	dloge = (math.log10(emax)-math.log10(emin))/(steps-1)
	for i in range(steps):
		energy = math.pow(10.0, math.log10(emin) + i*dloge)
		e_max.MeV(energy)
		flux_model = model.flux(e_min, e_max)
		flux_ff    = file_function.flux(e_min, e_max)
		if (abs(flux_model-flux_ff) > eps):
			print e_min, e_max, flux_model, flux_ff, flux_model-flux_ff
			n_violate += 1

	# Ramp up emax from above first node
	emin = 13.0
	e_min.MeV(emin)
	dloge = (math.log10(emax)-math.log10(emin))/(steps-1)
	for i in range(steps):
		energy = math.pow(10.0, math.log10(emin) + i*dloge)
		e_max.MeV(energy)
		flux_model = model.flux(e_min, e_max)
		flux_ff    = file_function.flux(e_min, e_max)
		if (abs(flux_model-flux_ff) > eps):
			print e_min, e_max, flux_model, flux_ff, flux_model-flux_ff
			n_violate += 1

	# Ramp down emin from above last node
	e_max.MeV(emax)
	dloge = (math.log10(emax)-math.log10(emin))/(steps-1)
	for i in range(steps):
		energy = math.pow(10.0, math.log10(emax) - i*dloge)
		e_min.MeV(energy)
		flux_model = model.flux(e_min, e_max)
		flux_ff    = file_function.flux(e_min, e_max)
		if (abs(flux_model-flux_ff) > eps):
			print e_min, e_max, flux_model, flux_ff, flux_model-flux_ff
			n_violate += 1

	# Ramp down emin from below last node
	emax = 873.0
	e_max.MeV(emax)
	dloge = (math.log10(emax)-math.log10(emin))/(steps-1)
	for i in range(steps):
		energy = math.pow(10.0, math.log10(emax) - i*dloge)
		e_min.MeV(energy)
		flux_model = model.flux(e_min, e_max)
		flux_ff    = file_function.flux(e_min, e_max)
		if (abs(flux_model-flux_ff) > eps):
			print e_min, e_max, flux_model, flux_ff, flux_model-flux_ff
			n_violate += 1

	# Print summary
	if n_violate == 0:
		print "compare_photon_flux: All tests satisfied accuracy."
	else:
		print "compare_photon_flux: "+str(n_violate)+" tests exceeded accuracy."

	# Return
	return


# ========================== #
# Compare energy flux values #
# ========================== #
def compare_energy_flux(model, file_function, emin=8.7, emax=1123.0, steps=13):
	"""
	Compare energy fluxes
	"""
	# Set accuracy
	eps       = 1.0e-5
	n_violate = 0
	
	# Allocate boundaries
	e_min = GEnergy()
	e_max = GEnergy()
	
	# Ramp up emax from below first node
	e_min.MeV(emin)
	dloge = (math.log10(emax)-math.log10(emin))/(steps-1)
	for i in range(steps):
		energy = math.pow(10.0, math.log10(emin) + i*dloge)
		e_max.MeV(energy)
		flux_model = model.eflux(e_min, e_max)
		flux_ff    = file_function.eflux(e_min, e_max)
		if (abs(flux_model-flux_ff) > eps):
			print e_min, e_max, flux_model, flux_ff, flux_model-flux_ff
			n_violate += 1

	# Ramp up emax from above first node
	emin = 13.0
	e_min.MeV(emin)
	dloge = (math.log10(emax)-math.log10(emin))/(steps-1)
	for i in range(steps):
		energy = math.pow(10.0, math.log10(emin) + i*dloge)
		e_max.MeV(energy)
		flux_model = model.eflux(e_min, e_max)
		flux_ff    = file_function.eflux(e_min, e_max)
		if (abs(flux_model-flux_ff) > eps):
			print e_min, e_max, flux_model, flux_ff, flux_model-flux_ff
			n_violate += 1

	# Ramp down emin from above last node
	e_max.MeV(emax)
	dloge = (math.log10(emax)-math.log10(emin))/(steps-1)
	for i in range(steps):
		energy = math.pow(10.0, math.log10(emax) - i*dloge)
		e_min.MeV(energy)
		flux_model = model.eflux(e_min, e_max)
		flux_ff    = file_function.eflux(e_min, e_max)
		if (abs(flux_model-flux_ff) > eps):
			print e_min, e_max, flux_model, flux_ff, flux_model-flux_ff
			n_violate += 1

	# Ramp down emin from below last node
	emax = 873.0
	e_max.MeV(emax)
	dloge = (math.log10(emax)-math.log10(emin))/(steps-1)
	for i in range(steps):
		energy = math.pow(10.0, math.log10(emax) - i*dloge)
		e_min.MeV(energy)
		flux_model = model.eflux(e_min, e_max)
		flux_ff    = file_function.eflux(e_min, e_max)
		if (abs(flux_model-flux_ff) > eps):
			print e_min, e_max, flux_model, flux_ff, flux_model-flux_ff
			n_violate += 1

	# Print summary
	if n_violate == 0:
		print "compare_energy_flux: All tests satisfied accuracy."
	else:
		print "compare_energy_flux: "+str(n_violate)+" tests exceeded accuracy."

	# Return
	return


# ======= #
# Test MC #
# ======= #
def test_mc(file_function, trials=10):
	"""
	Print a bunch of Monte Carlo energies.
	"""
	# Allocate boundaries
	e_min = GEnergy()
	e_max = GEnergy()
	ran   = GRan()

	# Set energy range
	e_min.TeV(0.1)
	e_max.TeV(100.0)
    
    # Loop over trials
	for i in range(trials):
		print file_function.mc(e_min, e_max, ran)

    # Return
    return
    

# ================ #
# Main entry point #
# ================ #
if __name__ == '__main__':
	"""
	Main entry point
	"""
	# Allocate power law
	model = GModelSpectralPlaw(5.7, -2.48)
	model["Prefactor"].scale(1.0e-16)
	model["PivotEnergy"].value(0.3)
	model["PivotEnergy"].scale(1.0e6)
	
	# Make file function
	file_function = make_file_function(model)
	
	# Compare photon fluxes
	compare_photon_flux(model, file_function)

	# Compare energy fluxes
	compare_energy_flux(model, file_function)
	
	# Test MC
	test_mc(file_function)
	