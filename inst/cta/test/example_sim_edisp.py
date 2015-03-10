#! /usr/bin/env python
# ==========================================================================
# This script simulates the energy dispersion distribution.
#
# If matplotlib is installed, the energy dispersion will be displayed on the
# screen.
#
# --------------------------------------------------------------------------
#
# Copyright (C) 2015 Florent Forest
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
import math


# ============== #
# Simulate Edisp #
# ============== #
def sim_edisp(edisp, etrue, eobs_max=40.0, ebins=500, nmc=100000):
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
        logEtrue = math.log10(etrue)

        # Create eobs axis
        eobs_axis = []
        dEobs    = eobs_max/ebins
        for i in range(ebins):
            eobs_axis.append((i+0.5)*dEobs)

        # Reset histogram
        counts = [0.0 for i in range(ebins)]

        # Allocate random number generator
        ran = gammalib.GRan()

        # Simulate observed energies
        for i in range(nmc):
            eobs     = edisp.mc(ran, logEtrue).TeV()
            index    = int(eobs/dEobs)
            if (index < ebins):
                counts[index] += 1.0

        # Create error bars
        error = [math.sqrt(c) for c in counts]

        # Get expected energy dispersion
        sum       = 0.0
        exp_edisp = []
        cumul     = []
        for i in range(ebins):
            eobs     = eobs_axis[i]
            if eobs <= edisp.ebounds_obs(logEtrue).emax().TeV():
                dLogEobs = math.log10(eobs_axis[i]) - math.log10(eobs_axis[i-1])
                value = edisp(math.log10(eobs), logEtrue) * dLogEobs * nmc
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
        #ax1.errorbar(eobs_axis, counts, error, ecolor='g')

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
    edisp = gammalib.GCTAEdispRmf("./caldb/dc1/rmf.fits")
    #edisp = gammalib.GCTAEdispPerfTable("./caldb/cta_dummy_irf.dat")
    #edisp = gammalib.GCTAEdisp2D("./caldb/edisp_matrix.fits")

    # Simulate Edisp
    sim_edisp(edisp, 20.0) # Energy in TeV
    #sim_edisp(edisp, 1.0) # Energy in TeV
