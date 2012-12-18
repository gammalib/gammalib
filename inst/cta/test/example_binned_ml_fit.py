#! /usr/bin/env python
# ==========================================================================
# This script performs a binned maximum likelihood analysis of CTA data.
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
from math import *


# =================== #
# CTA binned analysis #
# =================== #
def binned_analysis(model, cntmap, irf, caldb):
    """
    CTA binned analysis.
    """
    # Allocate empty observation container
    obs = GObservations()

    # Allocate empty CTA observation
    cta_obs = GCTAObservation()

    # Load counts map into CTA observation
    cta_obs.load_binned(cntmap)

    # Specify response for CTA observation
    cta_obs.response(irf, caldb)

    # Append CTA observation to observation container
    obs.append(cta_obs)

    # Load model to describe the data from XML file
    obs.models(model)

    # Allocate Levenberg-Marquardt optimizer
    opt = GOptimizerLM()

    # Optimize model parameters
    obs.optimize(opt)

    # Get maximum likelihood value
    logL = -(opt.value())

    # Get a copy of the model fitting results. We want a copy
    # here as the models are part of the observation container
    # "obs", and the container goes out of scope once the function
    # is left (and thus the models would also get out of scope)
    models = obs.models().copy()

    # Print optimizer results
    print(opt)

    # Create now a copy of the source model without background
    # only
    background = GModels()
    for m in models:
        if m.name() == "Background":
            background.append(m)

    # Assign background model for fitting
    obs.models(background)

    # Optimize background parameters
    obs.optimize(opt)

    # Get maximum likelihood value of background
    logL0 = -(opt.value())

    # Compute TS
    ts = 2.0 * (logL - logL0)

    # Print optimizer results
    print(opt)

    # Print TS
    print(" Test statistics ...........:", ts)

    # Return models
    return models


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    Example illustrating a binned data analysis.
    """
    # Dump header
    print("")
    print("*********************************************************")
    print("* Perform binned maximum likelihood fitting of CTA data *")
    print("*********************************************************")
    print("... please wait for a few seconds")

    # Set parameters
    irf    = "kb_E_50h_v3"
    caldb  = "../caldb"
    model  = "data/crab.xml"
    cntmap = "data/crab_cntmap.fits"

    # Perform binned analysis
    result = binned_analysis(model, cntmap, irf, caldb)

    # Print model results
    print(result)
