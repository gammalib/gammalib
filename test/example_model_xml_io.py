#!/usr/bin/env python
# ===========================================================================================#
# This script builds a set of models, saves the models into an XML file and loads the
# models back into memory.
# ===========================================================================================#

"""IO example:  GModels <--> XML file"""

from gammalib import *

#
# Set the location / centre of the source model to
# Right Ascension = 42.0 deg and Declination = 43 deg
# ===================================================
center = GSkyDir()
center.radec_deg(42.0, 43.0)

#
# Create a point source with power law spectral shape
# Normalization:   1.0e-7 ph/cm2/s/MeV @ 100 MeV
# Spectral index: -2.1
# ===================================================
point_spatial = GModelSpatialPtsrc(center)
point_spectrum = GModelSpectralPlaw(1.0e-7, -2.1)
point = GModelPointSource(point_spatial, point_spectrum)
point.name('My point source')

#
# Create a 2D Gaussian source with power law spectral shape
# Gaussian sigma:  3.0 deg
# Normalization:   4.2e-7 ph/cm2/s/MeV @ 100 MeV
# Spectral index: -2.4
# ===================================================
gauss_spatial = GModelRadialGauss(center, 3.0)
gauss_spectrum = GModelSpectralPlaw(4.2e-7, -2.4)
gauss = GModelExtendedSource(gauss_spatial, gauss_spectrum)
gauss.name('My Gaussian source')

#
# Allocate a GModels model container and append both source
# models to the model container. Print the content of the
# container
# =========================================================
models = GModels()
models.append(point)
models.append(gauss)
print models

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
models2 = GModels(filename)
print
print models2
