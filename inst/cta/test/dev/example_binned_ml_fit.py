#! /usr/bin/env python
# ==========================================================================
# This script performs a binned maximum likelihood analysis of CTA data.
#
# Copyright (C) 2012-2014 Juergen Knoedlseder
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
import time


# =================== #
# CTA binned analysis #
# =================== #
def binned_analysis(model, cntmap, irf, caldb):
    """
    Perform binned maximum likelihood fitting of CTA data using GammaLib
    classes. The function performs two fits, one with the full model and
    one with the background model only, to evaluate the Test Statistics
    value of the source component.
    """
    # Dump header
    print("")
    print("+===============================================+")
    print("| Binned maximum likelihood fitting of CTA data |")
    print("+===============================================+")

    # Get start CPU time
    tstart = time.clock()

    # Allocate empty observation container
    obs = gammalib.GObservations()

    # Allocate empty CTA observation
    cta_obs = gammalib.GCTAObservation()

    # Load counts map into CTA observation
    cta_obs.load(cntmap)

    # Specify response for CTA observation
    cta_obs.response(irf, caldb)

    # Append CTA observation to observation container
    obs.append(cta_obs)

    # Load model to describe the data from XML file
    obs.models(model)

    # Allocate Levenberg-Marquardt optimizer
    opt = gammalib.GOptimizerLM()

    # Optimize model parameters and compute errors
    obs.optimize(opt)
    obs.errors(opt)

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
    background = gammalib.GModels()
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
    print(" Test statistics ...........: %.3f" % ts)

    # Get stop CPU time
    tstop    = time.clock()
    telapsed = tstop - tstart
    print(" Elapsed time ..............: %.3f sec" % telapsed)

    # Return models
    return models


# ==================== #
# CTA stacked analysis #
# ==================== #
def stacked_analysis(model, cntmap, expcube, psfcube, bkgcube):
    """
    Perform stacked maximum likelihood fitting of CTA data using GammaLib
    classes. The function performs two fits, one with the full model and
    one with the background model only, to evaluate the Test Statistics
    value of the source component.
    """
    # Dump header
    print("")
    print("+================================================+")
    print("| Stacked maximum likelihood fitting of CTA data |")
    print("+================================================+")

    # Get start CPU time
    tstart = time.clock()

    # Allocate empty observation container
    obs = gammalib.GObservations()

    # Allocate empty CTA observation
    cta_obs = gammalib.GCTAObservation()

    # Load counts map into CTA observation
    cta_obs.load(cntmap)

    # Specify response for CTA observation
    exposure   = gammalib.GCTACubeExposure(expcube)
    psf        = gammalib.GCTACubePsf(psfcube)
    background = gammalib.GCTACubeBackground(bkgcube)
    cta_obs.response(exposure, psf, background)

    # Append CTA observation to observation container
    obs.append(cta_obs)

    # Load model to describe the data from XML file
    obs.models(model)

    # Allocate Levenberg-Marquardt optimizer
    opt = gammalib.GOptimizerLM()

    # Optimize model parameters and compute errors
    obs.optimize(opt)
    obs.errors(opt)

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
    background = gammalib.GModels()
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
    print(" Test statistics ...........: %.3f" % ts)

    # Get stop CPU time
    tstop    = time.clock()
    telapsed = tstop - tstart
    print(" Elapsed time ..............: %.3f sec" % telapsed)

    # Return models
    return models


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    Example illustrating binned data analyses.
    """
    # Dump header
    print("")
    print("**************************************************")
    print("* Perform maximum likelihood fitting of CTA data *")
    print("**************************************************")
    print("... please wait for a few seconds")

    # Set parameters
    irf     = "cta_dummy_irf"
    caldb   = "../caldb"
    model   = "data/crab.xml"
    cntmap  = "data/crab_cntmap.fits"
    expcube = "data/expcube.fits"
    psfcube = "data/psfcube.fits"
    bkgcube = "data/bkgcube.fits"

    # Perform binned analysis
    results_binned = binned_analysis(model, cntmap, irf, gammalib.GCaldb(caldb))

    # Perform stacked analysis
    results_cube = stacked_analysis(model, cntmap, expcube, psfcube, bkgcube)

    # Print model results
    #print(result)
