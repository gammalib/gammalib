#! /usr/bin/env python
# =====================================================================
# This script makes an image from a radial model and verifies that the
# model is correctly normalized.
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


# ============================ #
# Radial model testing routine #
# ============================ #
def test_radial_model(name, model, clobber=True):
    """
    Make an image of a spatial model and check that it integrates to 1.
    """
    # Make a 2D FITS image (cartesian projection, celestial coordinates)
    ra, dec, binsz, npix = 0, 0, 0.01, 300
    image = gammalib.GSkyMap("CAR", "CEL", ra, dec, -binsz, binsz, npix, npix, 1)

    # Fill the image
    energy = gammalib.GEnergy()
    time   = gammalib.GTime()
    for inx in range(image.npix()):
        dir        = image.inx2dir(inx)
        theta      = centre.dist(dir)
        image[inx] = model.eval(theta, energy, time)

    # Write it to file
    filename = name + '.fits'
    print('Writing', filename)
    image.save(filename, clobber)

    # Check that it integrates to 1
    integral = 0.0
    for inx in range(image.npix()):
        integral += image.solidangle(inx) * image[inx]
    print('integral = %g (should be 1)' % integral)
    print('integral error = %g (should be 0)' % (integral - 1))
    print('')


# ================ #
# Main entry point #
# ================ #
if __name__ == '__main__':
    """
    Test radial model integration.
    """
    # Set model location / centre
    centre = gammalib.GSkyDir()
    centre.radec_deg(0.3, 0.1)

    # Test models
    model = gammalib.GModelSpatialRadialGauss(centre, 0.3)
    test_radial_model(name='gauss', model=model)

    model = gammalib.GModelSpatialRadialDisk(centre, 0.8)
    test_radial_model(name='disk',  model=model)

    model = gammalib.GModelSpatialRadialShell(centre, 0.5, 0.1)
    test_radial_model(name='shell', model=model)
