#! /usr/bin/env python
# ==========================================================================
# This script tests the simulation of diffuse map cubes.
#
# If matplotlib is installed, the simulation results will be displayed.
#
# --------------------------------------------------------------------------
#
# Copyright (C) 2013 Juergen Knoedlseder
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


# ======================== #
# Simulate CTA observation #
# ======================== #
def simulate(xmlname, e_min, e_max, area, duration, dir, radius):
    """
    Simulate CTA observation.
    """
    # Allocate MC parameters
    emin = GEnergy()
    emax = GEnergy()
    tmin = GTime(0.0)
    tmax = GTime(duration)

    # Define MC parameters
    emin.TeV(e_min)
    emax.TeV(e_max)

    # Allocate random number generator
    ran = GRan()

    # Load models and extract first model
    models = GModels(xmlname)
    model = models[0]
    #print(model)
    #print(model.spatial())
    #print(model.spatial().spectrum())

    # Simulate photons
    photons = model.mc(area, dir, radius, emin, emax, tmin, tmax, ran)

    # Print photons
    print(str(len(photons)) + " photons simulated.")

    # Return photons
    return photons


# ============ #
# Show photons #
# ============ #
def show_photons(photons, xmlname, e_min, e_max, area, duration, dir, radius, ebins=30):
    """
    Show photons using matplotlib (if available).
    """
    # Only proceed if matplotlib is available
    try:
        # Import matplotlib
        import matplotlib.pyplot as plt

        # Create figure
        plt.figure(1)
        plt.title("MC simulated photon spectrum (" + str(e_min) + '-' + str(e_max) + " TeV)")

        # Setup energy range covered by data
        ebds = GEbounds(ebins, GEnergy(e_min, "TeV"), GEnergy(e_max, "TeV"))

        # Create energy axis
        energy = []
        for i in range(ebds.size()):
            energy.append(ebds.elogmean(i).TeV())

        # Fill histogram
        counts = [0.0 for i in range(ebds.size())]
        for photon in photons:
            index = ebds.index(photon.energy())
            counts[index] = counts[index] + 1.0

        # Create error bars
        error = [sqrt(c) for c in counts]

        # Get model values
        models = GModels(xmlname)
        mod    = models[0]
        model  = []
        t      = GTime()
        mod.spatial().set_mc_cone(dir, radius)
        for i in range(ebds.size()):
            eng    = ebds.elogmean(i)
            ewidth = ebds.ewidth(i)
            f      = mod.spatial().spectrum().eval(eng, t) * \
                     area * duration * ewidth.MeV()
            model.append(f)

        # Plot data
        plt.loglog(energy, counts, 'ro')
        plt.errorbar(energy, counts, error, fmt=None, ecolor='r')

        # Plot model
        plt.plot(energy, model, 'b-')

        # Set axes
        plt.xlabel("Energy (TeV)")
        plt.ylabel("Number of incident photons")

        # Create figure
        plt.figure(2)
        plt.title("MC simulated photon map")
        
        # Create RA and DEC arrays
        ra  = []
        dec = []
        for photon in photons:
            ra.append(photon.dir().ra_deg())
            dec.append(photon.dir().dec_deg())

        # Make scatter plot
        plt.scatter(ra, dec, marker=".")

        # Set axes
        plt.xlabel("Right Ascension (deg)")
        plt.ylabel("Declination (deg)")

        # Notify
        print("PLEASE CLOSE WINDOW TO CONTINUE ...")

        # Allocate histogram
        # Show plot
        plt.show()

    except ImportError:
        print("Matplotlib is not (correctly) installed on your system.")

    # Return
    return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    Simulate photons.
    """
    # Dump header
    print("")
    print("********************")
    print("* Simulate photons *")
    print("********************")

    # Set model filename
    xmlnames = ["data/model_diffuse_cube.xml"]

    # Set simulation parameters
    e_min    = 0.1                # 0.1 TeV
    e_max    = 100.0              # 100 TeV
    area     = 1.0                # 1 cm2
    duration = 1.0                # 1 sec
    dir      = GSkyDir()
    dir.radec_deg(84.17263, 22.01444)
    radius   = 1.0


    # Loop over models
    for xmlname in xmlnames:

        # Perform simulation
        photons = simulate(xmlname, e_min, e_max, area, duration, dir, radius)

        # Show photons
        show_photons(photons, xmlname, e_min, e_max, area, duration, dir, radius)
