#! /usr/bin/env python
# ==========================================================================
# This script tests the spatial map model GModelSpatialDiffuseMap.
#
# Copyright (C) 2012-2016 Juergen Knoedlseder
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


# ============== #
# Evaluate model #
# ============== #
def eval_model(model, binsz=0.02, nbins=200, ntrials=1000):
    """
    This function uses the eval method of the GModelSpatialDiffuseMap class to
    interpolate the skymap into a finer grid.
    """
    # Define centre and radius
    centre = gammalib.GSkyDir()
    centre.radec_deg(201.3651, -43.0191)
    radius = 1.5

    # Create CAR skymap
    map = gammalib.GSkyMap("CAR", "CEL", 201.3651, -43.0191, binsz, binsz, nbins, nbins)

    # Fill map
    sum = 0.0
    for i in range(map.npix()):
        dir       = map.inx2dir(i)
        if centre.dist_deg(dir) <= radius:
            photon    = gammalib.GPhoton(dir, gammalib.GEnergy(), gammalib.GTime())
            intensity = model.eval(photon) * map.solidangle(i)
            map[i]    = intensity
            sum      += intensity

    # Normalize map to number of trials
    norm = model.mc_norm(centre, radius)
    print(norm)
    print(sum)
    print(norm/sum)
    for i in range(map.npix()):
        map[i] = -map[i] * float(ntrials)/norm

    # Add MC simulated sky directions
    ran = gammalib.GRan()
    for i in range(ntrials):
        dir = model.mc(gammalib.GEnergy(), gammalib.GTime(), ran)
        inx = map.dir2inx(dir)
        if inx >= 0 and inx < map.npix():
            map[inx] += 1.0
        else:
            print('Pixel %d out of map boundary' % inx)

    # Determine map residuals
    residual = 0.0
    for i in range(map.npix()):
        residual += map[i]
    print('Residual: %e (%8.5f)' % (residual, residual/float(ntrials)))

    # Save map
    map.save('test_model_spatial_map.fits', True)

    # Return
    return


# ================ #
# Main entry point #
# ================ #
if __name__ == '__main__':

    # Allocate spatial map model
    model = gammalib.GModelSpatialDiffuseMap('data/cena_lobes_parkes.fits')
    print(model)

    # Evaluate skymap
    #eval_model(model)
    #eval_model(model, ntrials=10000000)
    eval_model(model, ntrials=10000000)
