#! /usr/bin/env python
# ==========================================================================
# Test energy integration of GCTABackground3D
#
# Copyright (C) 2017 Juergen Knoedlseder
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
import os
import math
import gammalib
import matplotlib.pyplot as plt


# ========================================= #
# Integrate using Trapezoid rule in log-log #
# ========================================= #
def integrate_trapezoid_loglog(bkg, logE1, logE2):
    """
    Integrate using Trapezoid rule in log-log
    """
    # Get function values
    logf1 = math.log(bkg(logE1, 0.0, 0.0))
    logf2 = math.log(bkg(logE2, 0.0, 0.0))
    
    # Compute ln of energy boundaries
    lnE1 = logE1 * gammalib.ln10
    lnE2 = logE2 * gammalib.ln10
    
    # Compute integral using Trapezoid rule in log-log
    integral = (lnE2 - lnE1) * (math.exp(lnE1 + logf1) +
                                math.exp(lnE2 + logf2)) * 0.5

    # Return integral
    return integral


# ========================================= #
# Integrate using Trapezoid rule in lin-lin #
# ========================================= #
def integrate_trapezoid_linlin(bkg, logE1, logE2):
    """
    Integrate using Trapezoid rule in lin-lin
    """
    # Set number of energy bins
    nbins = 1000

    # Compute log energy bin size
    dlogE  = (logE2-logE1)/float(nbins-1)

    # Perform integration
    integral = 0.0
    for i in range(nbins-1):
        logE_low  = logE1    + i*dlogE
        logE_high = logE_low + dlogE
        dE        = math.pow(10.0, logE_high) - math.pow(10.0, logE_low)
        integral += dE * (bkg(logE_low, 0.0, 0.0) + bkg(logE_high, 0.0, 0.0)) * 0.5

    # Return integral
    return integral


# =============== #
# Plot background #
# =============== #
def plot_background(bkg, emin=0.020, emax=120.0):
    """
    Plot background
    """
    # Generate energy vector
    energies = gammalib.GEnergies(100, gammalib.GEnergy(emin,'TeV'),
                                       gammalib.GEnergy(emax,'TeV'))

    # Set energy values
    x = [energy.TeV() for energy in energies]
    y = [bkg(energy.log10TeV(), 0.0, 0.0) for energy in energies]

    # Get node energies
    i_eng   = bkg.table().axis('ENERG')
    e_nodes = bkg.table().axis_nodes(i_eng)
    xo      = [math.pow(10.0,logE) for logE in e_nodes]
    yo      = [bkg(logE, 0.0, 0.0) for logE in e_nodes]

    # Compute integrals
    for i in range(len(xo)-1):
        if yo[i] > 0 and yo[i+1] > 0:
            integral_loglog = integrate_trapezoid_loglog(bkg, e_nodes[i], e_nodes[i+1])
            integral_linlin = integrate_trapezoid_linlin(bkg, e_nodes[i], e_nodes[i+1])
            print('%f-%f: %e %e (%f)' % (xo[i], xo[i+1], integral_loglog,
                                        integral_linlin,
                                        integral_loglog/integral_linlin))

    # Create figure
    plt.figure(1)
    plt.title('CTA Background')

    # Plot background
    plt.loglog(x, y, 'r-')
    plt.loglog(xo, yo, 'ro')

    # Set axes
    plt.xlabel("Energy (TeV)")
    plt.ylabel("Background (events/s/MeV/sr)")

    # Show plot
    plt.show()


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Load background
    bkg = gammalib.GCTABackground3D('../../caldb/data/cta/prod2/bcf/South_50h/irf_file.fits')
    print(bkg)

    # Plot background
    plot_background(bkg)
