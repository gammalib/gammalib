#! /usr/bin/env python
# ==========================================================================
# This script tests the accuracy of the GSkymap flux computation.
#
# --------------------------------------------------------------------------
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
import math
import gammalib

# Set integration paramters
dr = 0.001
da = 0.1
nr = 1000
na = 3600
dr_rad = dr * gammalib.deg2rad
da_rad = da * gammalib.deg2rad

# Allocate GSkymap
map = gammalib.GSkymap("CAR", "GAL", 0.0, 0.0, -1.0, 1.0, 360, 180)
#map = gammalib.GSkymap("GAL", 64, "RING")

# Loop over all pixels
for index in range(map.npix()):

    # Set map pixel flux
    map[index] = 1000000.0

    # Get flux using flux() method
    flux = map.flux(index)

    # Get sky direction of pixel
    dir = map.inx2dir(index)
    #print(dir)

    # Compute flux by integration
    sum = 0.0
    for ir in range(1,nr):
        r      = ir*dr
        r_rad  = r * gammalib.deg2rad
        sum_az = 0.0
        for ia in range(na):
            a   = ia*da
            loc = dir.copy()
            loc.rotate_deg(a,r)
            if map.dir2inx(loc) == index:
                sum_az += map(loc) * da_rad
        sum += sum_az * math.sin(r_rad) * dr_rad

    # Unset map pixel flux
    map[index] = 0.0

    # Print result
    print("%15s flux()=%.5f integration=%.5f delta=%.2f%%" % \
          (dir, flux, sum, (flux-sum)/sum*100.0))
