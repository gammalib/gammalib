#! /usr/bin/env python
# ===========================================================================================#
# This script illustrates how the GammaLib photon simulator works.
#
# Based on the MAGIC spectrum of the Crab nebula, and by assuming a powerlaw, it will
# create a Monte Carlo sample of photons.
#
# If matplotlib is installed, the spectrum will be displayed on the screen.
# ===========================================================================================#
from gammalib import *
from math import *


# ======================== #
# Simulate CTA observation #
# ======================== #
def simulate(xmlname, e_min, e_max, area, duration):
    """
    Simulate CTA observation.
    """
    # Allocate MC parameters
    dir = GSkyDir()
    emin = GEnergy()
    emax = GEnergy()
    tmin = GTime()
    tmax = GTime()

    # Define MC parameters
    dir.radec_deg(83.6331, 22.0145)
    radius = 10.0
    emin.TeV(e_min)
    emax.TeV(e_max)
    tmin.met(0.0)
    tmax.met(duration)

    # Allocate random number generator
    ran = GRan()

    # Load models and extract first model
    models = GModels(xmlname)
    model = models[0]
    print model

    # Simulate photons
    photons = model.mc(area, dir, radius, emin, emax, tmin, tmax, ran)

    # Print photons
    # for photon in photons:
    #	print photon
    print str(len(photons)) + " photons simulated."

    # Return photons
    return photons


# ============ #
# Show photons #
# ============ #
def show_photons(photons, xmlname, e_min, e_max, area, duration, ebins=30):
    """
    Show photons using matplotlib (if available).
    """
    # Only proceed if matplotlib is available
    try:
        # Import matplotlib
        import matplotlib.pyplot as plt

        # Create figure
        plt.figure(1)
        plt.title("MC simulated photon spectrum (" + str(e_min) + '-' + str(e_max) + " TeV)")

        # Setup energy range covered by data
        emin = GEnergy()
        emax = GEnergy()
        emin.TeV(e_min)
        emax.TeV(e_max)
        ebds = GEbounds()
        ebds.setlog(emin, emax, ebins)

        # Create energy axis
        energy = []
        for i in range(ebds.size()):
            energy.append(ebds.elogmean(i).TeV())

        # Fill histogram
        counts = [0.0 for i in range(ebds.size())]
        for photon in photons:
            index = ebds.index(photon.energy())
            counts[index] = counts[index] + 1.0

        # Create error bars
        error = [sqrt(c) for c in counts]

        # Get model values
        models = GModels(xmlname)
        crab = models[0]
        model = []
        d = GSkyDir()
        d.radec_deg(83.6331, 22.0145)
        t = GTime()
        for i in range(ebds.size()):
            eval = ebds.elogmean(i)
            ewidth = ebds.emax(i) - ebds.emin(i)
            f = crab.value(d, eval, t) * area * duration * ewidth.MeV()
            model.append(f)

        # Plot data
        plt.loglog(energy, counts, 'ro')
        plt.errorbar(energy, counts, error, fmt=None, ecolor='r')

        # Plot model
        plt.plot(energy, model, 'b-')

        # Set axes
        plt.xlabel("Energy (TeV)")
        plt.ylabel("Number of incident photons")

        # Notify
        print "PLEASE CLOSE WINDOW TO CONTINUE ..."

        # Allocate histogram
        # Show plot
        plt.show()

    except ImportError:
        print "Matplotlib is not (correctly) installed on your system."

    # Return
    return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    Simulate photons.
    """
    # Dump header
    print
    print "********************"
    print "* Simulate photons *"
    print "********************"

    # Set XML names
    xmlnames = ["data/crab.xml",
                "data/crab_file_function.xml",
                "data/crab_file_function_mod.xml"]

    # Set simulation parameters
    e_min = 0.1                # 0.1 TeV
    e_max = 100.0              # 100 TeV
    area = 3200000.0 * 1.0e4  # 3200000.0 m^2
    duration = 3600 * 5           # 5 hours

    # Loop over models
    for xmlname in xmlnames:

        # Perform simulation
        photons = simulate(xmlname, e_min, e_max, area, duration)

        # Show photons
        show_photons(photons, xmlname, e_min, e_max, area, duration)
