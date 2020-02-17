#! /usr/bin/env python
# ==========================================================================
# Test Npred computation for background models
#
# Copyright (C) 2020 Juergen Knoedlseder
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


# ===================== #
# Setup CTA observation #
# ===================== #
def setup_obs(dir, irf='South_50h', caldb='prod2', rad=2.0):
    """
    Setup observation
    """
    # Allocate CTA observation
    obs = gammalib.GCTAObservation()

    # Set CTA calibration database
    db = gammalib.GCaldb()
    if (gammalib.dir_exists(caldb)):
        db.rootdir(caldb)
    else:
        db.open('cta', caldb)

    # Set pointing direction for CTA observation
    pnt = gammalib.GCTAPointing()
    pnt.dir(dir)
    obs.pointing(pnt)

    # Set ROI
    roi    = gammalib.GCTARoi()
    centre = pnt.instdir(dir)
    roi.centre(centre)
    roi.radius(rad)

    # Set GTI
    ref = gammalib.GTimeReference(51544.5, 's', 'TT', 'LOCAL')
    gti = gammalib.GGti(ref)
    gti.append(gammalib.GTime(0.0, ref),
               gammalib.GTime(1800.0, ref))

    # Set energy boundaries
    ebounds = gammalib.GEbounds(gammalib.GEnergy(0.1, 'TeV'),
                                gammalib.GEnergy(100.0, 'TeV'))

    # Allocate event list
    events = gammalib.GCTAEventList()

    # Set ROI, GTI and energy boundaries for event list
    events.roi(roi)
    events.gti(gti)
    events.ebounds(ebounds)

    # Set the event list as the events for CTA observation
    obs.events(events)

    # Set instrument response for CTA observation
    obs.response(irf, db)

    # Set ontime, livetime, and deadtime correction factor for CTA observation
    obs.ontime(1800.0)
    obs.livetime(1800.0*0.98)
    obs.deadc(0.98)
    obs.id('000000')

    # Return CTA observation
    return obs


# ============================= #
# Setup radial acceptance model #
# ============================= #
def setup_model_radial():
    """
    """
    # Setup model
    pivot    = gammalib.GEnergy(1.0, 'TeV')
    spectral = gammalib.GModelSpectralPlaw(61.8e-6, -1.85, pivot)
    radial   = gammalib.GCTAModelRadialGauss(3.0)
    model    = gammalib.GCTAModelRadialAcceptance(radial, spectral)

    # Return model
    return model


# ====================== #
# Setup background model #
# ====================== #
def setup_model_background():
    """
    """
    # Setup model
    pivot    = gammalib.GEnergy(1.0, 'TeV')
    spectral = gammalib.GModelSpectralPlaw(61.8e-6, -1.85, pivot)
    spatial  = gammalib.GCTAModelRadialGauss(3.0)
    model    = gammalib.GCTAModelBackground(spatial, spectral)

    # Return model
    return model


# ================================== #
# Compute npred for background model #
# ================================== #
def npred(model, obs, dir):
    """
    Compute npred for background model
    """
    # Set energy and time
    energy = gammalib.GEnergy(1.0, "TeV")
    time   = gammalib.GTime()

    # Set RoI
    roi    = obs.roi()
    radius = roi.radius()
    centre = obs.pointing().instdir(dir)
    roi    = gammalib.GCTARoi(centre, radius)
    obs.events().roi(roi)

    # Compute npred
    npred = model.npred(energy, time, obs)

    # Return npred
    return npred


# ================================== #
# Compute npred for background model #
# ================================== #
def test():
    """
    """
    # Set pointing direction
    dir = gammalib.GSkyDir()
    dir.radec_deg(83.6331, 22.0145)
    
    # Get CTA observation
    obs = setup_obs(dir)

    # Get models
    radial     = setup_model_radial()
    background = setup_model_background()

    # Loop over RoI centre
    for i in range(10):

        # Set Right Ascension and Declination
        ra  = 83.6331 + i*0.1
        dec = 22.0145 + i*0.2
        dir.radec_deg(ra, dec)

        # Compute Npred
        value_radial     = npred(radial, obs, dir)
        value_background = npred(background, obs, dir)
        print(i, value_radial, value_background, (value_radial-value_background)/value_background)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Test
    test()
