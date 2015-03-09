#! /usr/bin/env python
# ==========================================================================
# This script simulates the energy dispersion distribution.
#
# If matplotlib is installed, the energy dispersion will be displayed on the
# screen.
#
# --------------------------------------------------------------------------
#
# Copyright (C) 2013-2015 Juergen Knoedlseder
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


# ============== #
# Simulate Edisp #
# ============== #
def sim_edisp(edisp, etrue, eobs_max=40.0, ebins=1000, nmc=100000):
    """
    Simulate Edisp and show results using matplotlib (if available).
    """
    # Check if matplotlib is available
    try:
        import matplotlib.pyplot as plt
        has_matplotlib = True
    except ImportError:
        print("Matplotlib is not (correctly) installed on your system.")
        has_matplotlib = False

    # Continue only if matplotlib is available
    if has_matplotlib:

        # Set log10(etrue)
        logEtrue = log10(etrue)

        # Create eobs axis
        eobs_axis = []
        dEobs    = eobs_max/ebins
        for i in range(ebins):
            eobs_axis.append((i+0.5)*dEobs)

        # Reset histogram
        counts = [0.0 for i in range(ebins)]

        # Allocate random number generator
        ran = GRan()

        # Simulate observed energies
        for i in range(nmc):
            eobs     = edisp.mc(ran, logEtrue).TeV()
            index    = int(eobs/dEobs)
            if (index < ebins):
                counts[index] += 1.0

        # Create error bars
        error = [sqrt(c) for c in counts]

        # Get expected energy dispersion
        sum = 0.0
        exp_edisp = []
        cumul = []
        for i in range(ebins):
            eobs     = eobs_axis[i]
            if eobs <= edisp.ebounds_obs(logEtrue).emax().TeV():
                dLogEobs = log10(eobs_axis[i]) - log10(eobs_axis[i-1])
                value = edisp(log10(eobs), logEtrue) * dLogEobs * nmc
            else:
                value = 0.0
            sum += value
            exp_edisp.append(value)
            cumul.append(sum/nmc)
        print(sum)

        # Create figure
        fig, ax1 = plt.subplots()

        # Plot simulated data
        ax1.plot(eobs_axis, counts, 'ro')
        ax1.errorbar(eobs_axis, counts, error, fmt=None, ecolor='g')

        # Plot energy dispersion
        ax1.plot(eobs_axis, exp_edisp, 'b-')

        # Set axes
        ax1.set_xlabel("Observed energy (TeV)")
        ax1.set_ylabel("Number of counts")

        ax2 = ax1.twinx()

        # Plot cumulative dispersion
        ax2.plot(eobs_axis, cumul, 'm')

        # Set axis
        ax2.set_ylabel("Cumulative probability")

        # Notify
        print("PLEASE CLOSE WINDOW TO CONTINUE ...")

        # Show plot
        plt.title("MC simulated Energy dispersion (true photon energy "+str(etrue)+" TeV)")
        plt.show()

    # Return
    return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    Simulate Energy dispersion.
    """
    # Dump header
    print("")
    print("******************************")
    print("* Simulate Energy dispersion *")
    print("******************************")

    # Load edisp
    edisp = GCTAEdispRmf("./caldb/dc1/rmf.fits")
    #edisp = GCTAEdispPerfTable("../caldb/cta_dummy_irf.dat")
    #edisp = GCTAEdisp2D("../caldb/data/cta/e/bcf/IFAE20120510_50h/irf_file_matrix.fits")

    # Simulate Edisp
    sim_edisp(edisp, 20.0)
