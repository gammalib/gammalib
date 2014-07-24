#! /usr/bin/env python
# ==========================================================================
# This script tests the model computations using the CTA response function.
#
# Copyright (C) 2011-2014 Jurgen Knodlseder
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


# ====================== #
# Set point source model #
# ====================== #
def ptsrc_model(ra=83.6331, dec=22.0145):
    """
    Set shell model.
    """
    # Set shell centre
    pos = gammalib.GSkyDir()
    pos.radec_deg(ra, dec)

    # Set spatial model
    spatial = gammalib.GModelSpatialPointSource(pos)

    # Set spectral model
    spectral = gammalib.GModelSpectralPlaw(1.0, -2.0, gammalib.GEnergy(100.0, "MeV"))

    # Set sky model
    model = gammalib.GModelSky(spatial, spectral)

    # Return model
    return model


# =============== #
# Set shell model #
# =============== #
def shell_model(ra=83.6331, dec=22.0145, radius=0.3, width=0.1):
    """
    Set shell model.
    """
    # Set shell centre
    centre = gammalib.GSkyDir()
    centre.radec_deg(ra, dec)

    # Set radial model
    radial = gammalib.GModelSpatialRadialShell(centre, radius, width, False)

    # Set spectral model
    spectral = gammalib.GModelSpectralPlaw(1.0, -2.0, gammalib.GEnergy(100.0, "MeV"))

    # Set sky model
    model = gammalib.GModelSky(radial, spectral)

    # Return model
    return model


# =============== #
# Set disk model #
# =============== #
def disk_model(ra=83.6331, dec=22.0145, radius=0.4):
    """
    Set disk model.
    """
    # Set disk centre
    centre = gammalib.GSkyDir()
    centre.radec_deg(ra, dec)

    # Set radial model
    radial = gammalib.GModelSpatialRadialDisk(centre, radius)

    # Set spectral model
    spectral = gammalib.GModelSpectralPlaw(1.0, -2.0, gammalib.GEnergy(100.0, "MeV"))

    # Set sky model
    model = gammalib.GModelSky(radial, spectral)

    # Return model
    return model


# ================== #
# Set Gaussian model #
# ================== #
def gauss_model(ra=83.6331, dec=22.0145, sigma=0.2):
    """
    Set Gaussian model.
    """
    # Set Gaussian centre
    centre = gammalib.GSkyDir()
    centre.radec_deg(ra, dec)

    # Set radial model
    radial = gammalib.GModelSpatialRadialGauss(centre, sigma)

    # Set spectral model
    spectral = gammalib.GModelSpectralPlaw(1.0, -2.0, gammalib.GEnergy(100.0, "MeV"))

    # Set sky model
    model = gammalib.GModelSky(radial, spectral)

    # Return model
    return model


# ========================= #
# Set Elliptical Disk model #
# ========================= #
def elliptical_model(ra=83.6331, dec=22.0145, semimajor=1.0, semiminor=0.5, posangle=45):
    """
    Set Elliptical Disk model.
    """
    # Set centre
    centre = gammalib.GSkyDir()
    centre.radec_deg(ra, dec)

    # Set elliptical model
    elliptical = gammalib.GModelSpatialEllipticalDisk(centre, semimajor, semiminor, posangle)

    # Set spectral model
    spectral = gammalib.GModelSpectralPlaw(1.0, -2.0, gammalib.GEnergy(100.0, "MeV"))

    # Set sky model
    model = gammalib.GModelSky(elliptical, spectral)

    # Return model
    return model


# ================= #
# Set diffuse model #
# ================= #
def diffuse_model(file="data/radio_map.fits"):
    """
    Set Diffuse model.
    """
    # Set spatial model
    spatial = gammalib.GModelSpatialDiffuseMap(file)

    # Set spectral model
    spectral = gammalib.GModelSpectralPlaw(1.0, -2.0, gammalib.GEnergy(100.0, "MeV"))

    # Set sky model
    model = gammalib.GModelSky(spatial, spectral)

    # Return model
    return model


# ========================== #
# Set binned CTA observation #
# ========================== #
def observation(ra=83.6331, dec=22.0145, emin=0.1, emax=100.0,
                binsz=0.02, npix=100, ebins=10,
                duration=1800.0, deadc=0.95):
    """
    Set binned CTA observation.
    """
    # Allocate observation
    obs = gammalib.GCTAObservation()

    # Set response
    #obs.response("cta_dummy_irf", GCaldb("../caldb"))
    exposure = gammalib.GCTAExposure("data/expcube.fits")
    psf      = gammalib.GCTAMeanPsf("data/psfcube.fits")
    obs.response(exposure, psf)

    # Set pointing
    dir = gammalib.GSkyDir()
    pnt = gammalib.GCTAPointing()
    dir.radec_deg(ra, dec)
    pnt.dir(dir)
    obs.pointing(pnt)

    # Set
    ebounds = gammalib.GEbounds(ebins, gammalib.GEnergy(emin, "TeV"), gammalib.GEnergy(emax, "TeV"))
    gti     = gammalib.GGti()
    tmin    = gammalib.GTime(0.0)
    tmax    = gammalib.GTime(duration)
    gti.append(tmin, tmax)
    map     = gammalib.GSkymap("CAR", "CEL", ra, dec, -binsz, binsz, npix, npix, ebins)
    cube    = gammalib.GCTAEventCube(map, ebounds, gti)
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
def test_irf(model, ra=83.6331, dec=22.0145, npix=100, filename="modmap.fits"):
    """
    Test IRF.
    """
    # Set CTA observation
    obs = observation(ra=ra, dec=dec, npix=npix)

    # Loop over all bins
    for bin in obs.events():
        irf = model.eval(bin, obs) * bin.size()
        bin.counts(irf)

    # Save observation
    obs.save(filename, True)

    # Print observation and model
    print(obs)
    print(model)

    # Return
    return


# ========================= #
# Test gradient computation #
# ========================= #
def test_grad(model, ra=83.6331, dec=22.0145, filename="gradmap.fits", ipar=2):
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
    Test model computation.
    """
    # Dump header
    print("")
    print("*************************")
    print("* Test model computation *")
    print("**************************")

    # Test IRF
    # test_irf(ptsrc_model())
    # test_irf(shell_model())
    # test_irf(disk_model())
    # test_irf(gauss_model())
    test_irf(elliptical_model())
    # test_irf(diffuse_model(), ra=201.3651, dec=-43.0191, npix=200)

    # Test Gradient
    # test_grad(ptsrc_model(), ipar=3)
    # test_grad(shell_model(), ipar=3)
    # test_grad(disk_model(), ipar=3)
    # test_grad(gauss_model(), ipar=3)
    # test_grad(diffuse_model(), ra=201.3651, dec=-43.0191, ipar=3)
