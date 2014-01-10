#! /usr/bin/env python
# ==========================================================================
# This script tests the Fermi-LAT analysis.
#
# Requires:
# - matplotlib
# - numpy
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

from gammalib import *
from math import *
import os


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
    lat.load_unbinned("data/p7v6/ft1.fits", "data/p7v6/ft2.fits", "data/p7v6/ltcube.fits")
    lat.response("P7SOURCE_V6", "../caldb")

    # Setup ROI covered by data
    instDir = GSkyDir()
    roi     = GLATRoi()
    instDir.radec_deg(83.6331, 22.0145)
    roi.centre(GLATInstDir(instDir))
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
    obs.optimize(opt)
    print(obs)


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
    lat.load_binned("data/p7v6/srcmap.fits", "data/p7v6/binned_expmap.fits", "data/p7v6/ltcube.fits")
    lat.response("P7SOURCE_V6", "../caldb")
    lat.response().force_mean(True)

    # Append LAT observation to container
    obs.append(lat)

    # Load models from XML file
    obs.models(xmlname)

    # Perform optimization
    log = GLog()
    log.cout(True)
    opt = GOptimizerLM(log)
    obs.optimize(opt)
    print(obs)
    print(opt)

    # Plot residuals
    try:
        plot_residuals(obs)
    except ImportError:
        print("Matplotlib is not (correctly) installed on your system.")


# ============= #
# Plot residual #
# ============= #
def plot_residuals(obs):
    """
    Plot residuals.
    """
    # Import
    import numpy as np
    import matplotlib.pyplot as plt

    # Set plot styles
    styles = ['b-', 'g-', 'y-', 'n-']

    # Open plot
    plt.figure(1)
    plt.title('Residuals')

    # Loop over observations
    for run in obs:

        # Get event cube
        cube = run.events()

        # Create energy axis
        energy = []
        ebds   = cube.ebounds()
        for i in range(ebds.size()):
            energy.append(ebds.elogmean(i).MeV())

        # Create spectrum
        counts = [0.0 for i in range(ebds.size())]
        for bin in cube:
            index = ebds.index(bin.energy())
            counts[index] = counts[index] + bin.counts()

        # Create error bars
        error = [sqrt(c) for c in counts]

        # Plot spectrum
        plt.loglog(energy, counts, 'ro', label='data')
        # plt.errorbar(energy, counts, error, fmt=None, ecolor='r')

        # Extract models
        sum_model = [0.0 for i in range(ebds.size())]
        for k, m in enumerate(obs.models()):
            model = [0.0 for i in range(ebds.size())]
            for i in range(cube.size()):
                bin = cube[i]
                index = ebds.index(bin.energy())
                prob = m.eval(bin, run)
                model[index] = model[index] + prob * bin.size()
            for i in range(ebds.size()):
                sum_model[i] = sum_model[i] + model[i]
            plt.loglog(energy, model, styles[k], label=m.name())
        plt.loglog(energy, sum_model, 'r-', label='total')

    # Put labels
    plt.xlabel("Energy (MeV)")
    plt.ylabel("Counts")
    plt.legend(loc="lower right")

    # Show residuals
    plt.show()

    # Return
    return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    Perform testing.
    """
    # Dump result
    print("")
    print("**************************************")
    print("* Fermi-LAT Python interface testing *")
    print("**************************************")

    # Analyse data
#    analyse_unbinned("data/p7v6/crab_model.xml")
    analyse_binned("data/p7v6/crab_model.xml")
