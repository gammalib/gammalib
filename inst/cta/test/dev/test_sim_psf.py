#! /usr/bin/env python
# ==========================================================================
# This script simulates the PSF distribution.
#
# If matplotlib is installed, the PSF will be displayed on the screen.
#
# --------------------------------------------------------------------------
#
# Copyright (C) 2013-2016 Juergen Knoedlseder
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
import math
import gammalib
try:
    import matplotlib.pyplot as plt
    has_matplotlib = True
except ImportError:
    print('Matplotlib is not (correctly) installed on your system.')
    has_matplotlib = False


# ============ #
# Simulate PSF #
# ============ #
def sim_psf(psf, energy, r_max=0.8, rbins=1000, nmc=1000000):
    """
    Simulate PSF and show results using matplotlib (if available).
    """
    # Set log10(energy)
    logE = math.log10(energy)

    # Create delta axis
    delta = []
    dr    = r_max/rbins
    for i in range(rbins):
        delta.append((i+0.5)*dr)

    # Reset histogram
    counts = [0.0 for i in range(rbins)]

    # Allocate random number generator
    ran = gammalib.GRan()

    # Simulate offsets
    for i in range(nmc):
        offset         = psf.mc(ran, logE) * gammalib.rad2deg
        index          = int(offset/dr)
        if (index < rbins):
            counts[index] += 1.0

    # Create error bars
    error = [math.sqrt(c) for c in counts]

    # Get expected PSF
    sum = 0.0
    values = []
    for i in range(rbins):
        r  = delta[i] * gammalib.deg2rad
        if r <= psf.delta_max(logE):
            value = psf(r, logE) * gammalib.twopi * math.sin(r) * \
                    dr * gammalib.deg2rad * nmc
        else:
            value = 0.0
        sum += value
        values.append(value)
    print(sum)

    # Continue only if matplotlib is available
    if has_matplotlib:

        # Create figure
        plt.figure(1)
        plt.title('MC simulated PSF ('+str(energy)+' TeV)')

        # Plot simulated data
        plt.plot(delta, counts, 'ro')
        plt.errorbar(delta, counts, error, ecolor='r')

        # Plot PSF
        plt.plot(delta, values, 'b-')
        
        # Set axes
        plt.xlabel('Offset angle (degrees)')
        plt.ylabel('Number of counts')

        # Notify
        print('PLEASE CLOSE WINDOW TO CONTINUE ...')

        # Show plot
        plt.show()

    # Return
    return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':

    # Dump header
    print('')
    print('****************')
    print('* Simulate PSF *')
    print('****************')

    # Load PSF
    #psf = gammalib.GCTAPsfKing('../caldb/data/cta/e/bcf/IFAE20120510_50h_King/irf_file.fits')
    #psf = gammalib.GCTAPsf2D('../caldb/data/cta/e/bcf/IFAE20120510_50h/irf_file.fits')
    #psf = gammalib.GCTAPsfPerfTable('../caldb/cta_dummy_irf.dat')
    psf = gammalib.GCTAPsfTable('../caldb/psf_table.fits[PSF_2D_TABLE]')

    # Simulate PSF
    sim_psf(psf, 1.0, r_max=0.5)
