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
import os



# ======================== #
# Simulate CTA observation #
# ======================== #
def simulate(observation,filename):
    """
    Simulate CTA observation.
    """

    # Allocate random number generator
    ran = GRan()

    # Build GCTAModelBackground with a diffuse cube    
    spat = GModelSpatialDiffuseCube(filename)
    pivot = GEnergy(1.0,"TeV")
    spec = GModelSpectralPlaw(1.0,0.0,pivot)
    bck = GCTAModelBackground(spat,spec)
    bck.instruments("CTA")

    # Simulate photons
    photons = bck.mc(observation, ran)

    # Print photons
    print(str(photons.size()) + " photons simulated.")
     
    # Return photons
    return photons,spat

# ============ #
# Show photons #
# ============ #
def show_photons(photons, mod,e_min,e_max,skydir,area,duration,ebins=30):#, e_min, e_max, area, duration, radius, ebins=30):
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
        #print counts
        #print error
        # Get model values
       
        model  = []
        t      = GTime()
        
        #update mc_cache  
        mod.set_mc_cone(skydir, radius)
        for i in range(ebds.size()):
            eng    = ebds.elogmean(i)
            ewidth = ebds.ewidth(i)
            f = mod.spectrum().eval(eng, t) * area * duration * ewidth.MeV()
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
        #plt.scatter(ra, dec, marker=".")
        x = skydir.ra_deg()
        y = skydir.dec_deg()
        
        xmin = x - sqrt(0.45) * radius
        xmax = x + sqrt(0.45) * radius
        ymin = y - sqrt(0.45) * radius
        ymax = y + sqrt(0.45) * radius
        
        plt.hist2d(ra,dec,bins=50,range=[[xmin,xmax],[ymin,ymax]])

        # Set axes
        plt.xlabel("Right Ascension (deg)")
        plt.ylabel("Declination (deg)")

        # Notify
        print("PLEASE CLOSE WINDOW TO CONTINUE ...")

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

    cube_filename = os.environ["GAMMALIB"]+"/test/data/test_cube.fits"
    duration = 180.0
    radius = 2.0
    skydir = GSkyDir()
    skydir.radec_deg(84.17263, 22.01444)
    e_min = 0.1
    e_max = 100.0

    # Create GCTAObservation  
    obs = GCTAObservation()
    
    # Create ROI
    roi = GCTARoi()
    instdir = GCTAInstDir()  
    instdir.dir(skydir)
    roi.centre(instdir)
    roi.radius(radius)
    
    # Create GTI
    tstart = GTime(0.)
    tstop = GTime(duration)
    gti = GGti()
    gti.append(tstart,tstop)
    
    # Create Ebounds
    emin = GEnergy(e_min,"TeV")
    emax = GEnergy(e_max,"TeV")
    ebds = GEbounds()
    ebds.append(emin,emax)
    
    #Create event list with Ebounds, ROI and GTI
    ev = GCTAEventList()
    ev.roi(roi)
    ev.gti(gti)    
    ev.ebounds(ebds)
    
    #Create Pointing
    pnt = GCTAPointing()
    pnt.dir(skydir)
    
    #Set observation parameters
    obs.ontime(gti.ontime())
    obs.livetime(gti.ontime())
    obs.deadc(1.0)
    obs.events(ev)
    obs.pointing(pnt)
        
    # simulate     
    photons,mod = simulate(obs,cube_filename)
    
    #Calculate surface area
    area = 4.0 * pi * (1.0-cos(radius/180.0*pi))
 
    # Show photons
    show_photons(photons, mod, e_min, e_max, skydir, area, duration)
