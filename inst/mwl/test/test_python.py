#! /usr/bin/env python

from gammalib import *
from math import *
import os
import numpy as np
import matplotlib.pyplot as plt


# ================================== #
# Extract data points from FITS file #
# ================================== #
def extract_data(filename):
    """
    Extract data from FITS file.
    """
    # Set plot styles
    ecolor = ['r',  'b',  'g']
    styles = ['ro', 'bo', 'go']

    # Open FITS file
    fits = GFits(filename)

    # Loop over all table HDUs
    k = 0
    for extno in range(fits.size()):
        if fits.hdu(extno).exttype() == GFitsHDU.HT_ASCII_TABLE or \
           fits.hdu(extno).exttype() == GFitsHDU.HT_BIN_TABLE:
            
            # Load spectrum
            spectrum = GMWLSpectrum()
            spectrum.load_fits(filename, extno)

            # Add spectral points
            x  = []
            y  = []
            dy = []
            for i in range(spectrum.size()):
                xval = spectrum.pointer(i).energy().MeV()
                yval = spectrum.pointer(i).flux()
                yerr = spectrum.pointer(i).flux_err()
                conv = xval*xval*1.6021765e-6
                x.append(xval)
                y.append(yval*conv)
                dy.append(yerr*conv)
            plt.loglog(x, y, styles[k])
            plt.errorbar(x, y, dy, fmt=None, ecolor=ecolor[k])
            
            # Increment style counter
            k = k + 1

    # Add axes
    plt.xlabel("Energy (MeV)")
    plt.ylabel("Flux (erg/cm2/s)")

    # Close FITS file
    fits.close()


# ============ #
# Fit spectrum #
# ============ #
def fit_spectrum(filename, xmlname):
    """
    Fit spectrum.
    """
    # Load observations
    obs     = GObservations()
    comptel = GMWLObservation(filename, "COMPTEL");
    lat     = GMWLObservation(filename, "LAT");
    hess    = GMWLObservation(filename, "HESS");
    obs.append(comptel)
    obs.append(lat)
    obs.append(hess)

    # Setup model
    models = GModels(xmlname)
    obs.models(models)

    # Perform optimization
    opt = GOptimizerLM()
    opt.max_iter(1000)
    obs.optimize(opt)
    print obs

    # Get model values
    x = np.power(10., np.arange(-1., 8., 0.1))
    y = []
    d = GSkyDir()
    e = GEnergy()
    t = GTime()
    for energy in x:
        e.MeV(energy)
        f = obs.models().value(d, e, t) * energy*energy*1.6021765e-6
        y.append(f)
    plt.plot(x, y)
    

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

    # Open plot
    plt.figure(1)
    plt.title('GammaLib fit of the Crab PWN')

    # Extract data from FITS file
    extract_data("data/crab_mwl.fits")

    # Fit spectrum
    fit_spectrum("data/crab_mwl.fits", "data/crab_mwl2.xml")

    # Show spectrum
    plt.show()
