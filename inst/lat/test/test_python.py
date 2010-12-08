#! /usr/bin/env python

from gammalib import *
from math import *
import os
import numpy as np
import matplotlib.pyplot as plt


# ===================== #
# Analyse unbinned data #
# ===================== #
def analyse_unbinned(xmlname):
    """
    Analyse unbinned LAT data.
    """
    # Allocate observations
    obs = GObservations()

    # Load LAT observation
    lat = GLATObservation()
    lat.load_unbinned("data/ft1.fits", "data/ft2.fits", "data/ltcube.fits");
    lat.response("Pass5_v0_Diffuse","irf");

    # Setup ROI covered by data
    instDir = GLATInstDir()
    roi     = GLATRoi()
    instDir.radec_deg(83.6331, 22.0145)
    roi.centre(instDir)
    roi.radius(7.5)
    lat.roi(roi)

    # Setup energy range covered by data
    emin = GEnergy()
    emax = GEnergy()
    emin.MeV(200.0)
    emax.MeV(20000.0)
    lat.ebounds().append(emin, emax)

    # Append LAT observation to container
    obs.append(lat)

    # Load models from XML file
    obs.models(xmlname)

    # Perform optimization
    opt = GOptimizerLM()
    opt.max_iter(1000)
    #obs.optimize(opt)
    print obs


# =================== #
# Analyse binned data #
# =================== #
def analyse_binned(xmlname):
    """
    Analyse binned LAT data.
    """
    # Allocate observations
    obs = GObservations()

    # Load LAT observation
    lat = GLATObservation()
    lat.load_binned("data/srcmap.fits", "data/binned_expmap.fits", "data/ltcube.fits");
    lat.response("Pass5_v0_Diffuse","irf");

    # Append LAT observation to container
    obs.append(lat)

    # Load models from XML file
    obs.models(xmlname)

    # Perform optimization
    log = GLog()
    log.cout(True)
    opt = GOptimizerLM(log)
    opt.max_iter(1000)
    obs.optimize(opt)
    print opt
    print obs
    
    # Plot residuals
    plot_residuals(obs)


# ============= #
# Plot residual #
# ============= #
def plot_residuals(obs):
    """
    Plot residuals.
    """
    # Set plot styles
    styles = ['b-', 'g-', 'y-', 'n-']

    # Open plot
    plt.figure(1)
    plt.title('Residuals')

    # Loop over observations
    for run in obs:
    
        # Get event cube
        cube = GLATEventCube(run.events())
    
        # Create energy axis
        energy = []
        ebds   = cube.ebds()
        for i in range(ebds.size()):
            energy.append(ebds.elogmean(i).MeV())
    
        # Create spectrum
        counts = [0.0 for i in range(ebds.size())]
        for bin in cube:
            index         = ebds.index(bin.energy())
            counts[index] = counts[index] + bin.counts()

        # Create error bars
        error = [sqrt(c) for c in counts]

        # Plot spectrum
        plt.loglog(energy, counts, 'ro')
        #plt.errorbar(energy, counts, error, fmt=None, ecolor='r')

        # Extract models
        sum_model = [0.0 for i in range(ebds.size())]
        for k, m in enumerate(obs.models()):
            model = [0.0 for i in range(ebds.size())]
            for i in range(cube.size()):
                bin          = cube.pointer(i)
                index        = ebds.index(bin.energy())
                prob         = m.eval_gradients(bin, run)
                model[index] = model[index] + prob * bin.size()
            for i in range(ebds.size()):
                sum_model[i] = sum_model[i] + model[i]
            plt.loglog(energy, model, styles[k])
        plt.loglog(energy, sum_model, 'r-')

    # Show residuals
    plt.show()


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    Perform testing.
    """
    # Dump result
    print
    print "****************************"
    print "* Python interface testing *"
    print "****************************"

    # Initialise success counter
    tests   = 0
    success = 0

    # Extract data from FITS file
    #extract_data("data/source_model.xml")

    # Analyse data
    #analyse_unbinned("data/source_model.xml")
    #analyse_binned("data/source_model.xml")   # Original
    #analyse_binned("data/source_model2.xml")  # Powerlaw for extragal. diffuse
    analyse_binned("data/source_model3.xml")  # No Crab, Powerlaw for extragal. diffuse

