#!/usr/bin/env python
# =====================================================================
# This script builds a set of models, saves the models into an XML file
# and loads the models back into memory.
#
# Copyright (C) 2011-2015 Juergen Knoedlseder
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
# =====================================================================
import gammalib

#
# Set the location / centre of the source model to
# Right Ascension = 42.0 deg and Declination = 43 deg
# ===================================================
centre = gammalib.GSkyDir()
centre.radec_deg(42.0, 43.0)

#
# Create a point source with power law spectral shape
# Normalization:   1.0e-7 ph/cm2/s/MeV @ 100 MeV
# Spectral index: -2.1
# ===================================================
point_spatial  = gammalib.GModelSpatialPointSource(centre)
point_spectrum = gammalib.GModelSpectralPlaw(1.0e-7, -2.1, gammalib.GEnergy(100.0, "MeV"))
point = gammalib.GModelSky(point_spatial, point_spectrum)
point.name('My point source')

#
# Create a 2D Gaussian source with power law spectral shape
# Gaussian sigma:  3.0 deg
# Normalization:   4.2e-7 ph/cm2/s/MeV @ 100 MeV
# Spectral index: -2.4
# ===================================================
gauss_spatial  = gammalib.GModelSpatialRadialGauss(centre, 3.0)
gauss_spectrum = gammalib.GModelSpectralPlaw(4.2e-7, -2.4, gammalib.GEnergy(100.0, "MeV"))
gauss = gammalib.GModelSky(gauss_spatial, gauss_spectrum)
gauss.name('My Gaussian source')

#
# Allocate a GModels model container and append both source
# models to the model container. Print the content of the
# container
# =========================================================
models = gammalib.GModels()
models.append(point)
models.append(gauss)
print(models)

#
# Write the model container to XML file and delete the container
# in-memory
# ==============================================================
filename = 'models.xml'
models.save(filename)
del models

#
# Now read the models back from the XML file into the container.
# Print again the content of the container.
# ==============================================================
models2 = gammalib.GModels(filename)
print("")
print(models2)
