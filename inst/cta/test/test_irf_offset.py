#! /usr/bin/env python
# ==========================================================================
# This script tests the offset angle dependence of the instrumental response
# function.
#
# It creates a set of models for different offset angles and evaluates the
# counts map for this set of models.
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
    spatial = GModelSpatialPtsrc(pos)

    # Set spectral model
    spectral = GModelSpectralPlaw(1.0, -2.0)

    # Set sky model
    model = GModelPointSource(spatial, spectral)

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
    centre = GSkyDir()
    centre.radec_deg(ra, dec)

    # Set radial model
    radial = GModelRadialShell(centre, radius, width, False)

    # Set spectral model
    spectral = GModelSpectralPlaw(1.0, -2.0)

    # Set sky model
    model = GModelExtendedSource(radial, spectral)

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
    centre = GSkyDir()
    centre.radec_deg(ra, dec)

    # Set radial model
    radial = GModelRadialDisk(centre, radius)

    # Set spectral model
    spectral = GModelSpectralPlaw(1.0, -2.0)

    # Set sky model
    model = GModelExtendedSource(radial, spectral)

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
    centre = GSkyDir()
    centre.radec_deg(ra, dec)

    # Set radial model
    radial = GModelRadialGauss(centre, sigma)

    # Set spectral model
    spectral = GModelSpectralPlaw(1.0, -2.0)

    # Set sky model
    model = GModelExtendedSource(radial, spectral)

    # Return model
    return model


# ========================== #
# Set binned CTA observation #
# ========================== #
def observation(ra=0.0, dec=0.0, binsz=0.05, npix=200, ebins=10):
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
    ebounds = GEbounds()
    emin    = GEnergy()
    emax    = GEnergy()
    emin.TeV(0.1)
    emax.TeV(100.0)
    ebounds.set_log(emin, emax, ebins)
    gti  = GGti()
    tmin = GTime(0.0)
    tmax = GTime(1800.0)
    gti.append(tmin, tmax)
    map = GSkymap("CAR", "CEL", ra, dec, -binsz, binsz, npix, npix, ebins)
    cube = GCTAEventCube(map, ebounds, gti)
    obs.events(cube)
    obs.ontime(1800.0)
    obs.livetime(1800.0*0.95)
    obs.deadc(0.95)

    # Optionally show observation
    # print(obs)

    # Return observation
    return obs


# ================ #
# Create model map #
# ================ #
def modmap(obs, models, phi=0, theta=0, filename="modmap.fits"):
    """
    Create model map by looping over all bins of the counts map and
    all models in the model container.
    """
    # Loop over all bins
    for bin in obs.events():

        # Compute IRF
        irf = 0.0
        for model in models:
            irf += model.eval(bin, obs) * bin.size()

        # Set bin
        bin.counts(irf)

    # Save observation
    obs.save(filename, True)

    # Return
    return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    Test offset angle dependence of IRF.
    """
    # Dump header
    print("")
    print("***************************************")
    print("* Test offset angle dependence of IRF *")
    print("***************************************")

    # Set set
    set = 2

    # Set CTA observation
    obs = observation()
    print obs

    # Set offset angle range
    # offsets = [0.0, 1.0, 2.0, 3.0]
    offsets = [0.0]

    # Loop over offset angles
    for offset in offsets:

        # Set models
        if set == 1:
            model1 = ptsrc_model(ra=0.0, dec=offset)
            model2 = ptsrc_model(ra=1.0, dec=0.0)
            model3 = ptsrc_model(ra=2.0, dec=0.0)
            model4 = ptsrc_model(ra=3.0, dec=0.0)
            model5 = ptsrc_model(ra=4.0, dec=0.0)
            models = [model1, model2, model3, model4, model5]
        elif set == 2:
            model1 = disk_model(ra=0.0, dec=offset)
            model2 = disk_model(ra=1.0, dec=0.0)
            model3 = disk_model(ra=2.0, dec=0.0)
            model4 = disk_model(ra=3.0, dec=0.0)
            model5 = disk_model(ra=4.0, dec=0.0)
            models = [model1, model2, model3, model4, model5]
        # model = shell_model(ra=0.0, dec=offset)
        # model = disk_model(ra=0.0, dec=offset)
        # model = gauss_model(ra=0.0, dec=offset)

        # Print model
        # print model

        # Set filename
        filename = "modmap_theta%2.2d.fits" % (int(offset * 10.0))

        # Create model map
        modmap(obs, models, phi=0.0, theta=0.0, filename=filename)
