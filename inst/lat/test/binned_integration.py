#! /usr/bin/env python
# ==========================================================================
# This script tests the binned integration.
#
# Requires:
# - matplotlib
# - numpy
#
# Copyright (C) 2012 Jurgen Knodlseder
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
        int = pow(100.0, gamma) / (1.0 - gamma) * \
            (pow(ebds.emax(i).MeV(), 1.0 - gamma) - pow(ebds.emin(i).MeV(), 1.0 - gamma))

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
            (pow(ebds.emin(i).MeV() / 100.0, -gamma) + pow(ebds.emax(i).MeV() / 100.0, -gamma))

        # Append integral to array
        integral.append(int)

    # Return array
    return integral


# ========================================= #
# Perform logarithmic Trapezoid integration #
# ========================================= #
def int_log_trapezoid(ebds, gamma):
    """
    Perform logarithmic Trapezoid integration.

    Parameters:
     ebds  - Energy boundaries
     gamma - Power law index
    """
    # Initialise result array
    integral = []

    # Loop over energy bins
    for i in range(ebds.size()):

        # Perform integration
        temp1 = ebds.emin(i).MeV() * pow(ebds.emin(i).MeV() / 100.0, -gamma)
        temp2 = ebds.emax(i).MeV() * pow(ebds.emax(i).MeV() / 100.0, -gamma)
        temp3 = log(ebds.emax(i).MeV() / ebds.emin(i).MeV())
        int_ = 0.5 * (temp1 + temp2) * temp3

        # Append integral to array
        integral.append(int_)

    # Return array
    return integral


# ========================================== #
# Perform arithmetic mean energy integration #
# ========================================== #
def int_mean(ebds, gamma):
    """
    Perform arithmetic mean energy integration.

    Parameters:
     ebds  - Energy boundaries
     gamma - Power law index
    """
    # Initialise result array
    integral = []

    # Loop over energy bins
    for i in range(ebds.size()):

        # Perform integration
        energy = 0.5 * (ebds.emax(i).MeV() + ebds.emin(i).MeV())
        int = pow(energy / 100.0, -gamma) * (ebds.emax(i).MeV() - ebds.emin(i).MeV())

        # Append integral to array
        integral.append(int)

    # Return array
    return integral


# =========================================== #
# Perform logaritmhic mean energy integration #
# =========================================== #
def int_logmean(ebds, gamma):
    """
    Perform logaritmhic mean energy integration.

    Parameters:
     ebds  - Energy boundaries
     gamma - Power law index
    """
    # Initialise result array
    integral = []

    # Loop over energy bins
    for i in range(ebds.size()):

        # Perform integration
        energy = ebds.elogmean(i).MeV()
        int = pow(energy / 100.0, -gamma) * (ebds.emax(i).MeV() - ebds.emin(i).MeV())

        # Append integral to array
        integral.append(int)

    # Return array
    return integral


# ========================================= #
# Perform geometric mean energy integration #
# ========================================= #
def int_geometric(ebds, gamma):
    """
    Perform geometric mean energy integration.

    Parameters:
     ebds  - Energy boundaries
     gamma - Power law index
    """
    # Initialise result array
    integral = []

    # Loop over energy bins
    for i in range(ebds.size()):

        # Perform integration
        energy = sqrt(ebds.emin(i).MeV() * ebds.emax(i).MeV())
        int = pow(energy / 100.0, -gamma) * (ebds.emax(i).MeV() - ebds.emin(i).MeV())

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
    gamma = 1.001
    ebins = 20

    # Setup energy boundaries
    ebds = set_ebds(ebins=ebins)

    # Create energy axis
    energy = []
    for i in range(ebds.size()):
        energy.append(ebds.elogmean(i).MeV())

    # Perform integrations
    analytical    = int_analytical(ebds, gamma)
    trapezoid     = int_trapezoid(ebds, gamma)
    log_trapezoid = int_log_trapezoid(ebds, gamma)
    mean          = int_mean(ebds, gamma)
    logmean       = int_logmean(ebds, gamma)
    geometric     = int_geometric(ebds, gamma)

    # Compute differences
    diff_trapezoid     = []
    diff_log_trapezoid = []
    diff_mean          = []
    diff_logmean       = []
    diff_geometric     = []
    for i in range(ebds.size()):
        diff_trapezoid.append(trapezoid[i] / analytical[i])
        diff_log_trapezoid.append(log_trapezoid[i] / analytical[i])
        diff_mean.append(mean[i] / analytical[i])
        diff_logmean.append(logmean[i] / analytical[i])
        diff_geometric.append(geometric[i] / analytical[i])

    # Open plot
    plt.figure(1)
    plt.title('Binned integration (Gamma=' + str(gamma) + ', ebins=' + str(ebins) + ', 0.1-20 GeV)')

    # Plot
    plt.semilogx(energy, diff_trapezoid, 'r-', label="Trapezoid")
    plt.semilogx(energy, diff_log_trapezoid, 'k-', label="Trapezoid (log)")
    plt.semilogx(energy, diff_mean, 'b-', label="Arithmetic mean")
    plt.semilogx(energy, diff_logmean, 'y-', label="Log mean")
    plt.semilogx(energy, diff_geometric, 'g.', label="Geometric mean")
    # plt.ylim(0.9,1.1)
    plt.xlabel('Energy (MeV)')
    plt.ylabel('Computed / analytical integral')
    plt.legend(loc="upper right")

    # Show plot
    plt.show()
