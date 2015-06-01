#! /usr/bin/env python
# ==========================================================================
# This script tests the accuracy of the HealPix flux computation. It sets
# a single pixel of the HealPix map, projects the HealPix map into a CAR
# map to compute the flux in the pixel in an indpendent way, and compare
# this flux estimates to that one derived from the flux() method.
#
# Note that this test relies on the correctness of the CAR flux()
# computation.
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
import gammalib

# Allocate HealPix map
nside = 16
npix  = 12 * nside*nside
print("Pixels %d" % npix)

# Loop over all pixels
for index in range(npix):

    # Allocate HealPix map
    map = gammalib.GSkymap("GAL", nside, "RING")

    # Create WCS map, set Healpix index and merge into WCS map
    wcs = gammalib.GSkymap("CAR", "GAL", 0.0, 0.0, -0.1, 0.1, 3600, 1800)
    map[index] = 1.0
    wcs += map

    # Get total flux in WCS map (in HealPix pixel only)
    total = 0.0
    for i in range(wcs.npix()):
        if map.dir2inx(wcs.inx2dir(i)) == index:
            total += wcs.flux(i)
        else:
            wcs[i] = 0

    # Save map
    wcs.save("healpix.fits", True)

    # Print results
    print("Index %d" % index)
    print("Flux in pixel from integration .....: %f" % total)
    print("Flux in pixel from flux() method ...: %f (%f%%)" % \
          (map.flux(index),(map.flux(index)-total)/total*100.0))