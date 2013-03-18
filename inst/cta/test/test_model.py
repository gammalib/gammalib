#! /usr/bin/env python
# ==========================================================================
# This script tests the model computations using the CTA response function.
#
# Copyright (C) 2011-2013 Jurgen Knodlseder
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


# ====================== #
# Set point source model #
# ====================== #
def ptsrc_model(ra=0.0, dec=0.0):
    """
    Set shell model.
    """
    # Set shell centre
    pos = GSkyDir()
    pos.radec_deg(ra, dec)

    # Set spatial model
    spatial = GModelSpatialPointSource(pos)

    # Set spectral model
    spectral = GModelSpectralPlaw(1.0, -2.0, 100.0)

    # Set sky model
    model = GModelSky(spatial, spectral)

    # Return model
    return model


# =============== #
# Set shell model #
# =============== #
def shell_model(ra=0.3, dec=0.3, radius=0.3, width=0.1):
    """
    Set shell model.
    """
    # Set shell centre
    center = GSkyDir()
    center.radec_deg(ra, dec)

    # Set radial model
    radial = GModelSpatialRadialShell(center, radius, width, False)

    # Set spectral model
    spectral = GModelSpectralPlaw(1.0, -2.0, 100.0)

    # Set sky model
    model = GModelSky(radial, spectral)

    # Return model
    return model


# =============== #
# Set disk model #
# =============== #
def disk_model(ra=359.6, dec=-0.2, radius=0.4):
    """
    Set disk model.
    """
    # Set disk centre
    center = GSkyDir()
    center.radec_deg(ra, dec)

    # Set radial model
    radial = GModelSpatialRadialDisk(center, radius)

    # Set spectral model
    spectral = GModelSpectralPlaw(1.0, -2.0, 100.0)

    # Set sky model
    model = GModelSky(radial, spectral)

    # Return model
    return model


# ================== #
# Set Gaussian model #
# ================== #
def gauss_model(ra=359.6, dec=+0.1, sigma=0.2):
    """
    Set Gaussian model.
    """
    # Set Gaussian centre
    center = GSkyDir()
    center.radec_deg(ra, dec)

    # Set radial model
    radial = GModelSpatialRadialGauss(center, sigma)

    # Set spectral model
    spectral = GModelSpectralPlaw(1.0, -2.0, 100.0)

    # Set sky model
    model = GModelSky(radial, spectral)

    # Return model
    return model


# ================= #
# Set diffuse model #
# ================= #
def diffuse_model(file="../../../test/data/cena_lobes_parkes.fits"):
    """
    Set Diffuse model.
    """
    # Set spatial model
    spatial = GModelSpatialDiffuseMap(file)

    # Set spectral model
    spectral = GModelSpectralPlaw(1.0, -2.0, 100.0)

    # Set sky model
    model = GModelSky(spatial, spectral)

    # Return model
    return model


# ========================== #
# Set binned CTA observation #
# ========================== #
def observation(ra=0.0, dec=0.0, emin=0.1, emax=100.0,
                binsz=0.02, npix=100, ebins=10,
                duration=1800.0, deadc=0.95):
    """
    Set binned CTA observation.
    """
    # Allocate observation
    obs = GCTAObservation()

    # Set response
    obs.response("kb_E_50h_v3", "../caldb")

    # Set pointing
    dir = GSkyDir()
    pnt = GCTAPointing()
    dir.radec_deg(ra, dec)
    pnt.dir(dir)
    obs.pointing(pnt)

    # Set
    ebounds = GEbounds(ebins, GEnergy(emin, "TeV"), GEnergy(emax, "TeV"))
    gti     = GGti()
    tmin    = GTime(0.0)
    tmax    = GTime(duration)
    gti.append(tmin, tmax)
    map     = GSkymap("CAR", "CEL", ra, dec, -binsz, binsz, npix, npix, ebins)
    cube    = GCTAEventCube(map, ebounds, gti)
    obs.events(cube)

    # Set ontime, livetime, and deadtime correction factor
    obs.ontime(duration)
    obs.livetime(duration * deadc)
    obs.deadc(deadc)

    # Optionally show observation
    # print(obs)

    # Return observation
    return obs


# ==================== #
# Test IRF computation #
# ==================== #
def test_irf(model, ra=0.0, dec=0.0, npix=100, filename="cntmap.fits"):
    """
    Test IRF.
    """
    # Print model
    print model

    # Set CTA observation
    obs = observation(ra=ra, dec=dec, npix=npix)

    # Loop over all bins
    for bin in obs.events():
        irf = model.eval(bin, obs) * bin.size()
        bin.counts(irf)

    # Save observation
    obs.save(filename, True)

    # Return
    return


# ========================= #
# Test gradient computation #
# ========================= #
def test_grad(model, ra=0.0, dec=0.0, filename="gradmap.fits", ipar=2):
    """
    Test gradient computation.
    """
    # Print model
    print(model)

    # Set CTA observation
    obs = observation(ra=ra, dec=dec, npix=npix)

    # Make sure that parameter is free
    model[ipar].free()

    # Loop over all bins
    for bin in obs.events():
        grad = obs.model_grad(model, bin, ipar)
        bin.counts(grad)

    # Save observation
    obs.save(filename, True)

    # Return
    return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    Test diffuse models.
    """
    # Dump header
    print("")
    print("***********************")
    print("* Test diffuse models *")
    print("***********************")

    # Test IRF
    # test_irf(ptsrc_model())
    # test_irf(shell_model())
    # test_irf(disk_model())
    # test_irf(gauss_model())
    test_irf(diffuse_model(), ra=201.3651, dec=-43.0191, npix=200)

    # Test Gradient
    # test_grad(ptsrc_model(), ipar=3)
    # test_grad(shell_model(), ipar=3)
    # test_grad(disk_model(), ipar=3)
    # test_grad(gauss_model(), ipar=3)
    # test_grad(diffuse_model(), ra=201.3651, dec=-43.0191, ipar=3)
