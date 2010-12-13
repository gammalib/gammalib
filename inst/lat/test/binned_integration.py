#! /usr/bin/env python

from gammalib import *
from math import *
import os
import numpy as np
import matplotlib.pyplot as plt


# ======================= #
# Setup energy boundaries #
# ======================= #
def set_ebds(ebins=20):
    """
    Setup energy boundaries.
    
    Keywords:
     ebins - Number of energy bins
    """
    # Setup energy range covered by data
    emin = GEnergy()
    emax = GEnergy()
    emin.MeV(100.0)
    emax.MeV(20000.0)
    ebds = GEbounds()
    ebds.setlog(emin, emax, ebins)

    # Return energy boundaries
    return ebds


# ============================== #
# Perform analytical integration #
# ============================== #
def int_analytical(ebds, gamma):
    """
    Perform analytical integration.
    
    Parameters:
     ebds  - Energy boundaries
     gamma - Power law index
    """
    # Initialise result array
    integral = []
    
    # Loop over energy bins
    for i in range(ebds.size()):

        # Perform analytical integration
        int = pow(100.0, gamma)/(1.0-gamma) * \
              (pow(ebds.emax(i).MeV(), 1.0-gamma) - pow(ebds.emin(i).MeV(), 1.0-gamma))
        
        # Append integral to array
        integral.append(int)
    
    # Return array
    return integral


# ============================= #
# Perform Trapezoid integration #
# ============================= #
def int_trapezoid(ebds, gamma):
    """
    Perform Trapezoid integration.
    
    Parameters:
     ebds  - Energy boundaries
     gamma - Power law index
    """
    # Initialise result array
    integral = []
    
    # Loop over energy bins
    for i in range(ebds.size()):

        # Perform integration
        int = 0.5 * (ebds.emax(i).MeV() - ebds.emin(i).MeV()) * \
              (pow(ebds.emin(i).MeV()/100.0, -gamma) + pow(ebds.emax(i).MeV()/100.0, -gamma))
        
        # Append integral to array
        integral.append(int)
    
    # Return array
    return integral


# =============================== #
# Perform Mean energy integration #
# =============================== #
def int_mean(ebds, gamma):
    """
    Perform mean energy integration.
    
    Parameters:
     ebds  - Energy boundaries
     gamma - Power law index
    """
    # Initialise result array
    integral = []

    # Loop over energy bins
    for i in range(ebds.size()):

        # Perform integration
        int = pow(ebds.elogmean(i).MeV()/100.0, -gamma) * \
              (ebds.emax(i).MeV() - ebds.emin(i).MeV())
        
        # Append integral to array
        integral.append(int)
    
    # Return array
    return integral


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    This script compares different energy integration methods used for binned
    analysis. The ScienceTools use the Trapezoid rule to integrate linearly in
    energy. The actual GammaLib implementation of the LAT interface use the
    logarithmic mean energy.
    """
    # Set parameters
    gamma = 2.0
    ebins = 20

    # Setup energy boundaries
    ebds = set_ebds(ebins=ebins)

    # Create energy axis
    energy = []
    for i in range(ebds.size()):
        energy.append(ebds.elogmean(i).MeV())

    # Perform integrations
    analytical = int_analytical(ebds, gamma)
    trapezoid  = int_trapezoid(ebds, gamma)
    mean       = int_mean(ebds, gamma)
    
    # Compute differences
    diff_trapezoid = []
    diff_mean      = []
    for i in range(ebds.size()):
        diff_trapezoid.append(trapezoid[i]/analytical[i])
        diff_mean.append(mean[i]/analytical[i])
    
    # Open plot
    plt.figure(1)
    plt.title('Binned integration')

    # Plot
    plt.semilogx(energy, diff_trapezoid, 'r-')
    plt.semilogx(energy, diff_mean, 'b-')
    #plt.ylim(0.9,1.1)
    plt.xlabel('Energy (MeV)')
    plt.ylabel('Computed / analytical integral')

    # Show plot
    plt.show()

