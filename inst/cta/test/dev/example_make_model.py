#! /usr/bin/env python
# ==========================================================================
# This script makes a model map from an XML file.
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


# ============================ #
# Radial model testing routine #
# ============================ #
def make_model_for_cntmap(cntmap, modname, xmlname, irf, caldb, clobber=True):
    """
    Make a model counts map from an XML model.
    """
    # Allocate empty CTA observation
    obs = gammalib.GCTAObservation()

    # Load counts map into CTA observation
    obs.load(cntmap)

    # Set response
    obs.response(irf, caldb)

    # Load models from XML file
    models = gammalib.GModels(xmlname)

    # Loop over all bins in counts map
    for event in obs.events():
        model = models.eval(event, obs) * event.size()
        event.counts(model)

    # Save CTA observation
    obs.save(modname, clobber)

    # Return
    return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    Example illustrating how to create a model map for a counts map.
    """
    # Dump header
    print("")
    print("***************************************************")
    print("* Create model map for counts map using XML model *")
    print("***************************************************")
    print("... please wait for a few seconds")

    # Set parameters
    irf     = "cta_dummy_irf"
    caldb   = "./caldb"
    xmlname = "data/crab.xml"
    cntmap  = "data/crab_cntmap.fits"
    modname = "model.fits"

    # Make the model
    make_model_for_cntmap(cntmap, modname, xmlname, irf, gammalib.GCaldb(caldb))
