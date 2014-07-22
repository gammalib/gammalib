#! /usr/bin/env python
# ==========================================================================
# This script performs a benchmark for maximum likelihood fitting of CTA
# data.
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


# ===================== #
# CTA unbinned analysis #
# ===================== #
def unbinned_analysis(model, evtfile, irf, caldb):
    """
    Perform unbinned maximum likelihood fitting of CTA data using GammaLib
    classes. The function performs two fits, one with the full model and
    one with the background model only, to evaluate the Test Statistics
    value of the source component.
    """
    # Dump header
    print("")
    print("+=================================================+")
    print("| Unbinned maximum likelihood fitting of CTA data |")
    print("+=================================================+")

    # Get start CPU time
    tstart = time.clock()

    # Allocate empty observation container
    obs = gammalib.GObservations()

    # Allocate empty CTA observation
    cta_obs = gammalib.GCTAObservation()

    # Load events into CTA observation
    cta_obs.load(evtfile)

    # Specify response for CTA observation
    cta_obs.response(irf, caldb)

    # Append CTA observation to observation container
    obs.append(cta_obs)

    # Load model to describe the data from XML file
    obs.models(model)

    # Allocate Levenberg-Marquardt optimizer
    opt = gammalib.GOptimizerLM()

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


# ======================= #
# CTA cube-style analysis #
# ======================= #
def cube_analysis(model, cntmap, expcube, psfcube):
    """
    Perform cube-style maximum likelihood fitting of CTA data using GammaLib
    classes. The function performs two fits, one with the full model and
    one with the background model only, to evaluate the Test Statistics
    value of the source component.
    """
    # Dump header
    print("")
    print("+===================================================+")
    print("| Cube-style maximum likelihood fitting of CTA data |")
    print("+===================================================+")

    # Get start CPU time
    tstart = time.clock()

    # Allocate empty observation container
    obs = gammalib.GObservations()

    # Allocate empty CTA observation
    cta_obs = gammalib.GCTAObservation()

    # Load counts map into CTA observation
    cta_obs.load(cntmap)

    # Specify response for CTA observation
    exposure = gammalib.GCTAExposure(expcube)
    psf      = gammalib.GCTAMeanPsf(psfcube)
    cta_obs.response(exposure, psf)

    # Append CTA observation to observation container
    obs.append(cta_obs)

    # Load model to describe the data from XML file
    obs.models(model)

    # Allocate Levenberg-Marquardt optimizer
    opt = gammalib.GOptimizerLM()

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
    Performs speed benchmark for CTA maximum likelihood fitting.
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
    evtfile = "data/crab_events.fits"
    cntmap  = "data/crab_cntmap.fits"
    expcube = "data/expcube.fits"
    psfcube = "data/psfcube.fits"

    # Perform unbinned analysis
    results_binned = unbinned_analysis(model, evtfile, irf, gammalib.GCaldb(caldb))

    # Perform binned analysis
    results_binned = binned_analysis(model, cntmap, irf, gammalib.GCaldb(caldb))

    # Perform cube-style analysis
    results_cube = cube_analysis(model, cntmap, expcube, psfcube)

    # Print model results
    #print(result)
