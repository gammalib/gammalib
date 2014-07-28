#! /usr/bin/env python
# ==========================================================================
# This script performs a benchmark of the response computation for all
# different analysis methods (unbinned, binned, stacked). Specifically it
# compares the total number of modelled counts to the expected number from
# Monte-Carlo simulations. This script is used to validate the absolute
# normalization of the CTA response computations.
#
# Copyright (C) 2014 Juergen Knoedlseder
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


# ============================ #
# CTA unbinned IRF computation #
# ============================ #
def unbinned_irf(model, events, irf, caldb, evtref):
    """
    Perform unbinned IRF computation.
    """
    # Dump header
    print(" Unbinned instrument response computation (Npred):")
    print(" =================================================")

    # Allocate empty CTA observation
    obs = gammalib.GCTAObservation()

    # Load events into CTA observation
    obs.load(events)

    # Specify response for CTA observation
    obs.response(irf, caldb)

    # Load model to describe the data from XML file
    models = gammalib.GModels(model)
    models.remove("Background")

    # Get start CPU time
    tstart = time.clock()

    # Compute Npred
    npred = obs.npred(models)

    # Get stop CPU time
    tstop      = time.clock()
    telapsed   = tstop - tstart
    difference = npred-evtref
    print("  Estimated number of counts : %.3f events" % npred)
    print("  Estimated - expected ......: %.3f events (%.1f%%)" % \
          (difference, difference/evtref*100.0))
    print("  Elapsed time ..............: %.3f sec" % telapsed)

    # Return
    return


# =========================== #
# CTA binned IRF computation #
# =========================== #
def binned_irf(model, cntmap, irf, caldb, cntref):
    """
    Perform binned IRF computation.
    """
    # Dump header
    print(" Binned instrument response computation:")
    print(" =======================================")

    # Allocate empty CTA observation
    obs = gammalib.GCTAObservation()

    # Load counts map into CTA observation
    obs.load(cntmap)

    # Specify response for CTA observation
    obs.response(irf, caldb)

    # Load model to describe the data from XML file
    models = gammalib.GModels(model)
    models.remove("Background")

    # Get start CPU time
    tstart = time.clock()

    # Compute the total content of the model cube
    total = 0.0
    for event in obs.events():
        total += models.eval(event, obs) * event.size()

    # Get stop CPU time
    tstop      = time.clock()
    telapsed   = tstop - tstart
    difference = total-cntref
    print("  Estimated number of counts : %.3f events" % total)
    print("  Estimated - expected ......: %.3f events (%.1f%%)" % \
          (difference, difference/cntref*100.0))
    print("  Elapsed time ..............: %.3f sec" % telapsed)

    # Return
    return


# =========================== #
# CTA stacked IRF computation #
# =========================== #
def stacked_irf(model, cntmap, expcube, psfcube, cntref):
    """
    Perform stacked IRF computation.
    """
    # Dump header
    print(" Stacked instrument response computation:")
    print(" ========================================")

    # Allocate empty CTA observation
    obs = gammalib.GCTAObservation()

    # Load counts map into CTA observation
    obs.load(cntmap)

    # Specify response for CTA observation
    exposure = gammalib.GCTAExposure(expcube)
    psf      = gammalib.GCTAMeanPsf(psfcube)
    obs.response(exposure, psf)

    # Load model to describe the data from XML file
    models = gammalib.GModels(model)
    models.remove("Background")

    # Get start CPU time
    tstart = time.clock()

    # Compute the total content of the model cube
    total = 0.0
    for event in obs.events():
        total += models.eval(event, obs) * event.size()

    # Get stop CPU time
    tstop      = time.clock()
    telapsed   = tstop - tstart
    difference = total-cntref
    print("  Estimated number of counts : %.3f events" % total)
    print("  Estimated - expected ......: %.3f events (%.1f%%)" % \
          (difference, difference/cntref*100.0))
    print("  Elapsed time ..............: %.3f sec" % telapsed)

    # Return
    return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    Perform CTA response computation benchmark.
    """
    # Dump header
    print("")
    print("**************************************")
    print("* CTA response computation benchmark *")
    print("**************************************")

    # Set response parameters
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

        # Perform unbinned computation
        unbinned_irf(model, events, irf, gammalib.GCaldb(caldb), evtref)

        # Perform binned computation
        binned_irf(model, cntmap, irf, gammalib.GCaldb(caldb), cntref)

        # Perform stacked computation
        stacked_irf(model, cntmap, expcube, psfcube, cntref)
