#! /usr/bin/env python
# ==========================================================================
# This script tests the GModelSpatialDiffuseConst model.
#
# Copyright (C) 2015 Juergen Knoedlseder
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
import gammalib
#import math


# ============== #
# Evaluate model #
# ============== #
def test_mc(samples=100000):
    """
    This function tests the mc() method.
    """
    # Create model
    model = gammalib.GModelSpatialDiffuseConst(1.0)

    # Create WCS map
    map = gammalib.GSkyMap("AIT", "GAL", 0.0, 0.0, -1.0, 1.0, 360, 180)

    # Set fixed parameters
    energy = gammalib.GEnergy()
    time   = gammalib.GTime()

    # Allocate random number generator
    ran = gammalib.GRan()

    # Fill map
    for i in range(samples):
        dir         = model.mc(energy, time, ran)
        pixel       = map.dir2inx(dir)
        map[pixel] += 1.0

    # Save map
    map.save("test_model_spatial_const.fits", True)

    # Return
    return


# ================ #
# Main entry point #
# ================ #
if __name__ == '__main__':
    """
    Main entry point
    """
    # Test Monte Carlo method
    test_mc()
