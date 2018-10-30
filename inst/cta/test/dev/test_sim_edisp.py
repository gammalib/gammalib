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
import time


# ============== #
# Simulate Edisp #
# ============== #
def sim_edisp(edisp, etrue, theta=0.0, ebins=1000, nmc=1000000):
    """
    Simulate Edisp and show results using matplotlib (if available).
    """
    # Check if matplotlib is available
    try:
        import matplotlib.pyplot as plt
        has_matplotlib = True
    except ImportError:
        print('Matplotlib is not (correctly) installed on your system.')
        has_matplotlib = False

    # Continue only if matplotlib is available
    if has_matplotlib:

        # Create ereco axis
        ereco_axis   = []
        ereco_bounds = edisp.ereco_bounds(etrue, theta)
        ereco_min    = ereco_bounds.emin().TeV()
        ereco_max    = ereco_bounds.emax().TeV()
        dE           = (ereco_max-ereco_min)/ebins
        for i in range(ebins):
            ereco_axis.append(ereco_min+i*dE)
        print('Ereco: %.3f - %.3f TeV' % (ereco_min, ereco_max))

        # Create etrue axis
        etrue_axis   = []
        etrue_bounds = edisp.etrue_bounds(etrue, theta)
        etrue_min    = etrue_bounds.emin().TeV()
        etrue_max    = etrue_bounds.emax().TeV()
        dEtrue       = (etrue_max-etrue_min)/ebins
        for i in range(ebins):
            etrue_axis.append(etrue_min+i*dEtrue)
        print('Etrue: %.3f - %.3f TeV' % (etrue_min, etrue_max))

        # Reset histogram
        counts = [0.0 for i in range(ebins)]

        # Allocate random number generator
        ran = gammalib.GRan()

        # Simulate observed energies
        tstart = time.clock()
        for i in range(nmc):
            eobs = edisp.mc(ran, etrue, theta).TeV()
            index = int((eobs-ereco_min)/dE-0.5)
            if (index < ebins):
                counts[index] += 1.0
        telapsed = (time.clock() - tstart)
        print('Elapsed CPU time (sec) %f' % telapsed)

        # Dump probabilities
        print('Probability integrated over Ereco: %f' % \
              edisp.prob_erecobin(ereco_bounds.emin(), ereco_bounds.emax(), etrue, theta))

        # Create error bars
        error = [math.sqrt(c) for c in counts]

        # Get expected energy dispersion
        edisp_ereco = []
        sum_ereco   = 0.0
        for i in range(ebins):
            ereco = gammalib.GEnergy(ereco_axis[i], 'TeV')
            if ereco.TeV() >= ereco_min and ereco.TeV() <= ereco_max:
                value = edisp(ereco, etrue, theta) * dE * 1.0e6
            else:
                value = 0.0
            sum_ereco += value
            edisp_ereco.append(value)
        print('Sum computed over Ereco: %f' % sum_ereco)

        # Get energy dispersion versus etrue
        edisp_etrue = []
        sum_etrue   = 0.0
        for i in range(ebins):
            ereco = gammalib.GEnergy(etrue_axis[i], 'TeV')
            if ereco.TeV() >= etrue_min and ereco.TeV() <= etrue_max:
                value = edisp(etrue, ereco, theta) * dEtrue * 1.0e6
            else:
                value = 0.0
            sum_etrue += value
            edisp_etrue.append(value)
        print('Sum computed over Etrue: %f' % sum_etrue)

        # Normalise expectation
        if sum_ereco > 0.0:
            norm = nmc/sum_ereco
            for i in range(ebins):
                edisp_ereco[i] *= norm

        # Create figure
        fig = plt.figure(figsize=(14,6))

        # Add title
        title = r'MC simulated energy dispersion (E = %s, theta = %.2f)' % \
                (str(etrue), theta*gammalib.rad2deg)
        fig.suptitle(title, fontsize=16)

        # Plot energy dispersion versus ereco
        ax1 = fig.add_subplot(121)
        ax1.plot(ereco_axis, edisp_ereco, '-', linewidth=2)
        ax1.step(ereco_axis, counts, 'r', where='mid', linewidth=1)
        ymin, ymax = ax1.get_ylim()
        ax1.plot([etrue.TeV(), etrue.TeV()], [ymin,ymax], 'g-', linewidth=2)
        ax1.set_xlabel('Reconstructed energy (TeV)')
        ax1.set_ylabel('Number of counts')
        ax1.set_title(r'Energy dispersion versus E$_{\rm reco}$')

        # Plot energy dispersion versus etrue
        ax2 = fig.add_subplot(122)
        ax2.plot(etrue_axis, edisp_etrue, '-', linewidth=2)
        ymin, ymax = ax2.get_ylim()
        ax2.plot([etrue.TeV(), etrue.TeV()], [ymin,ymax], 'g-', linewidth=2)
        ax2.set_xlabel('True energy (TeV)')
        ax2.set_ylabel('Arbitrary unit')
        ax2.set_title(r'Energy dispersion versus E$_{\rm true}$')

        # Notify
        print('PLEASE CLOSE WINDOW TO CONTINUE ...')

        # Show plot
        plt.show()

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Dump header
    print('')
    print('******************************')
    print('* Simulate Energy dispersion *')
    print('******************************')

    # Load edisp
    #edisp = gammalib.GCTAEdispRmf('./caldb/dc1/rmf.fits')
    #edisp = gammalib.GCTAEdispPerfTable('./caldb/cta_dummy_irf.dat')

    # Simulate 2D energy dispersion matrix
    #edisp = gammalib.GCTAEdisp2D('../caldb/edisp_matrix.fits')
    #sim_edisp(edisp, gammalib.GEnergy(0.1, 'TeV'))
    #sim_edisp(edisp, gammalib.GEnergy(1.0, 'TeV'))
    #sim_edisp(edisp, gammalib.GEnergy(10.0, 'TeV'))
    #sim_edisp(edisp, gammalib.GEnergy(1.0, 'TeV'), theta=5.0*gammalib.deg2rad)

    # Test H.E.S.S. energy dispersion
    edisp = gammalib.GCTAEdisp2D('/project-data/hess/hess_dl3_dr1/data/hess_dl3_dr1_obs_id_023523.fits')
    #sim_edisp(edisp, gammalib.GEnergy(0.2, 'TeV'))
    #sim_edisp(edisp, gammalib.GEnergy(0.3, 'TeV'))
    #sim_edisp(edisp, gammalib.GEnergy(0.4, 'TeV'), theta=2.0*gammalib.deg2rad)
    #sim_edisp(edisp, gammalib.GEnergy(0.5, 'TeV'), theta=2.5*gammalib.deg2rad)
    #sim_edisp(edisp, gammalib.GEnergy(0.7, 'TeV'))
    #sim_edisp(edisp, gammalib.GEnergy(1.0, 'TeV'))
    #sim_edisp(edisp, gammalib.GEnergy(0.7, 'TeV'))
    sim_edisp(edisp, gammalib.GEnergy(5.0, 'TeV'))
    #sim_edisp(edisp, gammalib.GEnergy(50.0, 'TeV'))
    #sim_edisp(edisp, gammalib.GEnergy(1.0, 'TeV'), theta=2.5*gammalib.deg2rad)


