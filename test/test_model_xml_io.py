#!/usr/bin/env python
"""IO example:  GModels <--> XML file"""

from gammalib import *

# Make a GModels object containing two models
center = GSkyDir()
center.radec_deg(42, 43)

point_spatial = GModelSpatialPtsrc(center)
point_spectrum = GModelSpectralPlaw(1.0e-7, -2.1)
point = GModelPointSource(point_spatial, point_spectrum)
point.name('point')

gauss_spatial = GModelRadialGauss(center, 3)
gauss_spectrum = GModelSpectralPlaw(4.2e-7, -2.4)
gauss = GModelExtendedSource(gauss_spatial, gauss_spectrum)
gauss.name('gauss')

models = GModels()
models.append(point)
models.append(gauss)
print models

# Write GModels to file and delete it in-memory
filename = 'models.xml'
models.save(filename)
del models

# Now read it back
models2 = GModels(filename)
print models2
