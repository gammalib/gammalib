#! /usr/bin/env python
# ===========================================================================================#
# This script displays the radial acceptance model that may be used to model the
# CTA radial acceptance.
#
# Required 3rd party modules:
# - matplotlib
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
# ===========================================================================================#
import matplotlib.pyplot as plt
from gammalib import *
from math import *


# ========== #
# Show model #
# ========== #
def show_model(xmlfile):
	"""
	Show radial acceptance model from the XML file using matplotlib.
	"""
	# Load the model
	models = GModels(xmlfile)
	
	# Extract radial acceptance model
	radial = cast_GCTAModelRadialAcceptance(models["Background"]).radial()
	
	# Create angular axis (from 0 to 4 deg)
	thetas = [i*0.05 for i in range(80)]
	
	# Get model values
	values      = [radial.eval(theta) for theta in thetas]
	values_grad = [radial.eval_gradients(theta) for theta in thetas]
	
	# Create figure
	plt.figure(1)
	plt.title("Radial acceptance model ("+radial.type()+")")
				
	# Plot data
	plt.plot(thetas, values,      'r-')
	plt.plot(thetas, values_grad, 'ro')
		
	# Set axes
	plt.xlabel("Offset angle (deg)")
	plt.ylabel("Function value")
		
	# Show plot
	plt.show()

	# Return
	return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
	"""
	Show radial acceptance models.
	"""
	# Dump header
	print
	print "*********************************"
	print "* Show radial acceptance models *"
	print "*********************************"

    # Display various models
	show_model("data/crab.xml")
	show_model("data/crab_poly.xml")
	show_model("data/crab_profile.xml")
