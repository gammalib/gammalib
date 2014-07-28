#! /usr/bin/env python
# ==========================================================================
# This script performs a likelihood profile computation
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
try:
    import matplotlib.pyplot as plt
    has_matplotlib = True
except ImportError:
    print("Matplotlib is not (correctly) installed on your system.")
    has_matplotlib = False


# ======================= #
# Show likelihood profile #
# ======================= #
def show_likelihood_profile(result):
    """
    """
    # Continue only if matplotlib is available
    if has_matplotlib:

        # Create log-likelohood profile
        plt.figure(1)
        plt.title("Log-likelihood profile for parameter \""+result["parname"]+"\"")

        # Plot logL profile
        plt.plot(result["value"], result["logL"], 'ro-')

        # Build gradient array
        gradient = []
        for value in result["value"]:
            grad = result["opt"]["gradient"] * (value - result["opt"]["value"]) + \
                   result["opt"]["logL"]
            gradient.append(grad)

        # Plot gradient
        plt.plot(result["value"], gradient, 'b-')

        # Set axes
        plt.xlabel(result["parname"])
        plt.ylabel("log-likelihood")

        # Create Npred profile
        plt.figure(2)
        plt.title("Npred profile for parameter \""+result["parname"]+"\"")

        # Plot Npred profile
        plt.plot(result["value"], result["npred"], 'ro-')

        # Set axes
        plt.xlabel(result["parname"])
        plt.ylabel("Npred")

        # Notify
        print("PLEASE CLOSE WINDOW TO CONTINUE ...")

        # Show plot
        plt.show()

    # Return
    return


# ========================== #
# Compute likelihood profile #
# ========================== #
def compute_likelihood_profile(obs, parname="RA", scale=1.0e-4, steps=10):
    """
    Compute the log-likelihood profile for a given parameter.
    
    Keywords:
     parname - Parameter name
     scale   - Variation size of parameter for profile computation
     steps   - Number of steps to the left and to the right of the optimum
    """
    # Free parameter for which profile is to be determined
    obs.models()[0][parname].free()
    obs.models()[0]["DEC"].free()

    # Allocate Levenberg-Marquardt optimizer
    log = gammalib.GLog()
    log.cout(True)
    opt = gammalib.GOptimizerLM(log)

    # Optimize model parameters
    obs.optimize(opt)

    # Log optimizer into console
    print(opt)
    print(obs)
    print(obs.models())

    # Get optimizer parameter and gradient
    opt_logL     = obs.logL()
    opt_npred    = obs.npred()
    opt_value    = obs.models()[0][parname].value()
    opt_gradient = obs.models()[0][parname].gradient()

    # Log results into console
    print("%.10f %.10f %.10f %.10f" % (opt_value, opt_logL, opt_npred, opt_gradient))

    # Allocate arrays
    values = []
    logLs  = []
    npreds = []

    # Loop over values
    for i in range(-steps, steps+1):
    
        # Set new model parameter
        obs.models()[0][parname].remove_range()
        obs.models()[0][parname].value(opt_value+i*scale)
        
        # Evaluate log-likelihood and retrieve results
        obs.eval()
        value = obs.models()[0][parname].value()
        logL  = obs.logL()
        npred = obs.npred()
        values.append(value)
        logLs.append(logL)
        npreds.append(npred)
        
        # Log results into console
        print("%.10f %.10f %.10f" % (value, obs.logL(), obs.npred()))

    # Build result dictionary
    result = {'parname': parname,
              'opt': {'logL': opt_logL, 'npred': opt_npred, \
                      'value': opt_value, 'gradient': opt_gradient}, \
              'value': values, 'logL': logLs, 'npred': npreds}

    # Return result
    return result


# ===================================================== #
# CTA unbinned or binned likelihood profile computation #
# ===================================================== #
def original_likelihood_profile(model, cntmap, irf, caldb, cntref):
    """
    Perform binned likelihood profile computation.
    """
    # Dump header
    print("Unbinned or Binned likelihood profile computation:")
    print("==================================================")

    # Allocate empty observation container
    obs = gammalib.GObservations()

    # Allocate empty CTA observation
    cta = gammalib.GCTAObservation()

    # Load counts map into CTA observation
    cta.load(cntmap)

    # Specify response for CTA observation
    cta.response(irf, caldb)

    # Append CTA observation to observation container
    obs.append(cta)

    # Load model
    models = gammalib.GModels(model)
    obs.models(models)

    # Get start CPU time
    tstart = time.clock()

    # Compute likelihood profile
    result = compute_likelihood_profile(obs)

    # Get stop CPU time
    tstop      = time.clock()
    telapsed   = tstop - tstart
    print(" Elapsed time ..............: %.3f sec" % telapsed)

    # Show likelihood profile
    show_likelihood_profile(result)

    # Return
    return


# ========================================== #
# CTA stacked likelihood profile computation #
# ========================================== #
def stacked_likelihood_profile(model, cntmap, expcube, psfcube, cntref):
    """
    Perform stacked likelihood profile computation.
    """
    # Dump header
    print("Stacked likelihood profile computation:")
    print("=======================================")

    # Allocate empty observation container
    obs = gammalib.GObservations()

    # Allocate empty CTA observation
    cta = gammalib.GCTAObservation(cntmap, expcube, psfcube)

    # Append CTA observation to observation container
    obs.append(cta)

    # Load model
    models = gammalib.GModels(model)
    obs.models(models)

    # Get start CPU time
    tstart = time.clock()

    # Compute likelihood profile
    result = compute_likelihood_profile(obs)

    # Get stop CPU time
    tstop      = time.clock()
    telapsed   = tstop - tstart
    print(" Elapsed time ..............: %.3f sec" % telapsed)

    # Show likelihood profile
    show_likelihood_profile(result)

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
    print("* CTA likelihood profile computation *")
    print("**************************************")

    # Set test
    #test = "point"
    test = "gauss"
    #test = "disk"
    #test = "ldisk"
    #test = "shell"
    #test = "ellipse"
    #test = "diffuse"

    # Set response parameters
    irf     = "cta_dummy_irf"
    caldb   = "../caldb"
    #irf     = "irf_file.fits"
    #caldb   = "../caldb/data/cta/e/bcf/IFAE20120510_50h"
    expcube = "data/expcube.fits"
    psfcube = "data/psfcube.fits"

    # Set test dependent filenames
    if test == "point":
        model  = "data/crab_ptsrc.xml"
        events = "data/crab_events.fits"
        cntmap = "data/crab_cntmap.fits"
        cntref = 934.3 # +/- 3.1 (from ctobssim)
    elif test == "disk":   
        model  = "data/crab_disk.xml"
        events = "data/crab_disk_events.fits"
        cntmap = "data/crab_disk_cntmap.fits"
        cntref = 933.6 # +/- 3.1 (from ctobssim)
    elif test == "ldisk":   
        model  = "crab_ldisk.xml"
        events = "crab_ldisk_events.fits"
        cntmap = "data/crab_disk_cntmap.fits"
        cntref = 933.6 # +/- 3.1 (from ctobssim)
    elif test == "gauss":   
        model  = "data/crab_gauss.xml"
        events = "data/crab_gauss_events.fits"
        cntmap = "data/crab_gauss_cntmap.fits"
        cntref = 935.1 # +/- 3.1 (from ctobssim)
    elif test == "shell":   
        model  = "data/crab_shell.xml"
        events = "data/crab_shell_events.fits"
        cntmap = "data/crab_shell_cntmap.fits"
        cntref = 933.5 # +/- 3.1 (from ctobssim)
    elif test == "ellipse":   
        model  = "data/crab_edisk.xml"
        events = "data/crab_edisk_events.fits"
        cntmap = "data/crab_edisk_cntmap.fits"
        cntref = 845.9 # +/- 2.9 (from ctobssim)
    elif test == "diffuse":
        model  = "data/radio.xml"
        events = "data/radio_events.fits"
        cntmap = "data/radio_cntmap.fits"
        cntref = 337.5 # +/- 1.8 (from ctobssim)

    # Perform unbinned computation
    original_likelihood_profile(model, events, irf, gammalib.GCaldb(caldb), cntref)

    # Perform binned computation
    #original_likelihood_profile(model, cntmap, irf, gammalib.GCaldb(caldb), cntref)

    # Perform stacked computation
    #stacked_likelihood_profile(model, cntmap, expcube, psfcube, cntref)
