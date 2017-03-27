#! /usr/bin/env python
# ==========================================================================
# Test spectral model
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
import math
import gammalib
import matplotlib.pyplot as plt
import numpy as np


# ============== #
# Evaluate model #
# ============== #
def eval_model(model, emin=0.030, emax=100.0, ntrials=1000):
    """
    Evaluates model.
    """
    # Set energy boundaries
    e_min = gammalib.GEnergy(emin, 'TeV')
    e_max = gammalib.GEnergy(emax, 'TeV')

    # Generate model spectrum
    sum     = 0.0
    flux    = model.flux(e_min, e_max)
    norm    = float(ntrials)/flux
    x       = []
    y       = []
    ebounds = gammalib.GEbounds(100, e_min, e_max)
    for i in range(len(ebounds)):
        energy = ebounds.elogmean(i)
        value  = model.eval(energy) * norm * ebounds.ewidth(i).MeV()
        sum   += value
        if value > 1.0:
            x.append(energy.TeV())
            y.append(value)
    print(sum, ntrials)
    
    # Generate Monte Carlo spectrum
    energies = []
    ran      = gammalib.GRan()
    time     = gammalib.GTime()
    for i in range(ntrials):
        energy = model.mc(e_min, e_max, time, ran)
        energies.append(energy.TeV())

    # Create figure
    plt.figure()

    # Create histogram
    plt.hist(energies, bins=np.logspace(math.log10(emin), math.log10(emax), 50),
             log=True)
    plt.gca().set_xscale("log")

    # Show model
    plt.plot(x, y, 'r-')

    # Set axes
    plt.xlabel('Energy (TeV)')
    plt.ylabel('Arbitrary units')

    # Show spectrum
    plt.show()

    # Return
    return


# ================ #
# Main entry point #
# ================ #
if __name__ == '__main__':

    # Allocate spectral model
    model = gammalib.GModelSpectralSmoothBrokenPlaw(5.7e-16,
                                                    -2.48,
                                                    gammalib.GEnergy(0.3, 'TeV'),
                                                    -3.5,
                                                    gammalib.GEnergy(1.0, 'TeV'),
                                                    0.2)

    # Evaluate model
    eval_model(model, ntrials=1000000)
