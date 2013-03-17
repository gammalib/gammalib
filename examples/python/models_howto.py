#! /usr/bin/env python
# ==========================================================================
# Scope
#
#   This script provides an example for creating and handling models with
#   GammaLib
#
# Usage
#   ./model_howto.py
#
# -------------------------------------------------------------------------
#
# Copyright (C) 2013 Juergen Knoedlseder
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


# =========================================== #
# Create a model container filled with models #
# =========================================== #
def create_models():
    """
    """
    # Create model container
    models = GModels()
    
    # Create a powerlaw model for the Crab
    crabdir = GSkyDir()
    crabdir.radec_deg(83.6331, 22.0145)
    spatial  = GModelSpatialPointSource(crabdir)
    spectral = GModelSpectralPlaw(5.7e-16, -2.48, 100.0)
    model    = GModelSky(spatial, spectral)
    models.append(model)
    
    # Return models
    return models


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    Example script for creating and handling model.
    """
    # Dump header
    print("")
    print("******************************")
    print("* Example for model handling *")
    print("******************************")
    print("")

    # Create XML document
    models = create_models()
    print models
    