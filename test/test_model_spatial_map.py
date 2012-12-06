#! /usr/bin/env python
# ==========================================================================
# This script tests the GModelSpatialMap spatial map model.
#
# Copyright (C) 2012 Juergen Knoedlseder
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


# ============== #
# Evaluate model #
# ============== #
def eval_model(model, binsz=0.01, nbins=1000):
    """
    This function uses the eval method of the GModelSpatialMap class to
    interpolate the skymap into a finer grid.
    """
    # Create CAR skymap
    map = GSkymap("CAR", "CEL", 201.3651, -43.0191, binsz, binsz, nbins, nbins)

    # Fill map
    for i in range(map.npix()):
        dir = map.pix2dir(i)
        intensity = model.eval(dir)
        map[i] = intensity

    # Save map
    map.save("test_model_spatial_map.fits")

    # Return
    return


# ================ #
# Main entry point #
# ================ #
if __name__ == '__main__':
    """
    Main entry point
    """
    # Allocate spatial map model
    model = GModelSpatialMap("data/cena_lobes_parkes.fits")
    print model

    # Evaluate skymap
    eval_model(model)
