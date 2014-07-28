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
    classes.
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

    # Print optimizer results
    print(opt)
    print(obs.models())

    # Get stop CPU time
    tstop    = time.clock()
    telapsed = tstop - tstart
    print(" Elapsed time ..............: %.3f sec" % telapsed)

    # Return
    return


# =================== #
# CTA binned analysis #
# =================== #
def binned_analysis(model, cntmap, irf, caldb):
    """
    Perform binned maximum likelihood fitting of CTA data using GammaLib
    classes.
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

    # Print optimizer results
    print(opt)
    print(obs.models())

    # Get stop CPU time
    tstop    = time.clock()
    telapsed = tstop - tstart
    print(" Elapsed time ..............: %.3f sec" % telapsed)

    # Return
    return


# ==================== #
# CTA stacked analysis #
# ==================== #
def stacked_analysis(model, cntmap, expcube, psfcube):
    """
    Perform stacked maximum likelihood fitting of CTA data using GammaLib
    classes.
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

    # Print optimizer results
    print(opt)
    print(obs.models())

    # Get stop CPU time
    tstop    = time.clock()
    telapsed = tstop - tstart
    print(" Elapsed time ..............: %.3f sec" % telapsed)

    # Return
    return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    Performs speed benchmark for CTA maximum likelihood fitting.
    """
    # Dump header
    print("")
    print("********************************************")
    print("* CTA maximum likelihood fitting benchmark *")
    print("********************************************")

    # Set parameters
    irf     = "cta_dummy_irf"
    caldb   = "../caldb"
    expcube = "data/expcube.fits"
    psfcube = "data/psfcube.fits"

    # Set tests
    tests = ["point", "disk", "gauss", "shell", "ellipse", "diffuse"]

    # Loop over tests
    for test in tests:

        # Set test dependent filenames
        if test == "point":
            model  = "data/crab_ptsrc.xml"
            events = "data/crab_events.fits"
            cntmap = "data/crab_cntmap.fits"
            cntref = 934.3 # +/- 3.1 (from ctobssim)
            evtref = 934.3 # +/- 3.1 (from ctobssim)
        elif test == "disk":   
            model  = "data/crab_disk.xml"
            events = "data/crab_disk_events.fits"
            cntmap = "data/crab_disk_cntmap.fits"
            cntref = 933.6 # +/- 3.1 (from ctobssim)
            evtref = 933.6 # +/- 3.1 (from ctobssim)
        elif test == "gauss":   
            model  = "data/crab_gauss.xml"
            events = "data/crab_gauss_events.fits"
            cntmap = "data/crab_gauss_cntmap.fits"
            cntref = 935.1 # +/- 3.1 (from ctobssim)
            evtref = 935.1 # +/- 3.1 (from ctobssim)
        elif test == "shell":   
            model  = "data/crab_shell.xml"
            events = "data/crab_shell_events.fits"
            cntmap = "data/crab_shell_cntmap.fits"
            cntref = 933.5 # +/- 3.1 (from ctobssim)
            evtref = 933.5 # +/- 3.1 (from ctobssim)
        elif test == "ellipse":   
            model  = "data/crab_edisk.xml"
            events = "data/crab_edisk_events.fits"
            cntmap = "data/crab_edisk_cntmap.fits"
            cntref = 845.9 # +/- 2.9 (from ctobssim)
            evtref = 845.9 # +/- 2.9 (from ctobssim)
        elif test == "diffuse":
            model  = "data/radio.xml"
            events = "data/radio_events.fits"
            cntmap = "data/radio_cntmap.fits"
            cntref = 337.5 # +/- 1.8 (from ctobssim)
            evtref = 368.0 # +/- 1.9 (from ctobssim)

        # Print header
        print("")
        print("Model: %s" % (test))

        # Perform unbinned analysis
        unbinned_analysis(model, events, irf, gammalib.GCaldb(caldb))

        # Perform binned analysis
        binned_analysis(model, cntmap, irf, gammalib.GCaldb(caldb))

        # Perform stacked analysis
        stacked_analysis(model, cntmap, expcube, psfcube)
