#! /usr/bin/env python
# ==========================================================================
# This script tests the event simulation using GCTAModelBackground.
#
# If matplotlib is installed, the simulation results will be displayed.
#
# --------------------------------------------------------------------------
#
# Copyright (C) 2013 Michael Mayer
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
def simulate(observation, filename):
    """
    Simulate CTA observation.
    """

    # Allocate random number generator
    ran = GRan()

    # Build GCTAModelBackground with a diffuse cube    
    spat  = GModelSpatialDiffuseCube(filename)
    pivot = GEnergy(1.0,"TeV")
    spec  = GModelSpectralPlaw(1.0,0.0,pivot)
    bck   = GCTAModelBackground(spat,spec)
    bck.instruments("CTA")

    # Simulate events
    events = bck.mc(observation, ran)

    # Print events
    print(str(events.size()) + " events simulated.")
     
    # Return events
    return events, spat


# =========== #
# Show events #
# =========== #
def show_events(events, mod, e_min, e_max, skydir, duration, ebins=30):
    """
    Show events using matplotlib (if available).
    """
    # Only proceed if matplotlib is available
    try:
        # Import matplotlib
        import matplotlib.pyplot as plt

        # Create figure
        plt.figure(1)
        plt.title("MC simulated event spectrum (" + str(e_min) + '-' + str(e_max) + " TeV)")
        
        # Setup energy range covered by data
        ebds = GEbounds(ebins, GEnergy(e_min, "TeV"), GEnergy(e_max, "TeV"))

        # Create energy axis
        energy = []
        for i in range(ebds.size()):
            energy.append(ebds.elogmean(i).TeV())

        # Fill histogram
        counts = [0.0 for i in range(ebds.size())]
        for event in events:
            index = ebds.index(event.energy())
            counts[index] = counts[index] + 1.0
        
        # Create error bars
        error = [sqrt(c) for c in counts]
        #print counts
        #print error

        # Get model values       
        model  = []
        t      = GTime()
        mod.set_mc_cone(skydir, radius)
        for i in range(ebds.size()):
            eng    = ebds.elogmean(i)
            ewidth = ebds.ewidth(i)
            value  = mod.spectrum().eval(eng, t) * duration * ewidth.MeV()
            model.append(value)

        # Plot data
        plt.loglog(energy, counts, 'ro')
        plt.errorbar(energy, counts, error, fmt=None, ecolor='r')

        # Plot model
        plt.plot(energy, model, 'b-')

        # Set axes
        plt.xlabel("Energy (TeV)")
        plt.ylabel("Number of events")

        # Create figure
        plt.figure(2)
        plt.title("MC simulated event map")
        
        # Create RA and DEC arrays
        ra  = []
        dec = []
        for event in events:
            ra.append(event.dir().dir().ra_deg())
            dec.append(event.dir().dir().dec_deg())

        # Make 2D histogram (or scatter plot in case that the 2D histogram
        # is not available)
        if hasattr(plt, 'hist2d'):
            x    = skydir.ra_deg()
            y    = skydir.dec_deg()
            xmin = x - sqrt(0.45) * radius
            xmax = x + sqrt(0.45) * radius
            ymin = y - sqrt(0.45) * radius
            ymax = y + sqrt(0.45) * radius
            plt.hist2d(ra,dec,bins=50,range=[[xmin,xmax],[ymin,ymax]])
        else:
            plt.scatter(ra, dec, marker=".")
 
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
    Simulate events.
    """
    # Dump header
    print("")
    print("*******************")
    print("* Simulate events *")
    print("*******************")

    #cube_filename = os.environ["GAMMALIB"]+"/test/data/test_cube.fits"
    cube_filename = "../../../test/data/test_cube.fits"
    duration      = 1.0
    radius        = 2.0
    e_min         = 0.1
    e_max         = 100.0
    skydir        = GSkyDir()
    skydir.radec_deg(84.17263, 22.01444)

    # Create GCTAObservation  
    obs = GCTAObservation()
    
    # Create ROI
    roi     = GCTARoi()
    instdir = GCTAInstDir()  
    instdir.dir(skydir)
    roi.centre(instdir)
    roi.radius(radius)
    
    # Create GTI
    tstart = GTime(0.)
    tstop  = GTime(duration)
    gti    = GGti()
    gti.append(tstart,tstop)
    
    # Create Ebounds
    emin = GEnergy(e_min,"TeV")
    emax = GEnergy(e_max,"TeV")
    ebds = GEbounds()
    ebds.append(emin,emax)
    
    # Create event list with Ebounds, ROI and GTI
    ev = GCTAEventList()
    ev.roi(roi)
    ev.gti(gti)    
    ev.ebounds(ebds)
    
    # Create Pointing
    pnt = GCTAPointing()
    pnt.dir(skydir)
    
    # Set observation parameters
    obs.ontime(gti.ontime())
    obs.livetime(gti.ontime())
    obs.deadc(1.0)
    obs.events(ev)
    obs.pointing(pnt)
        
    # Simulate events  
    events, mod = simulate(obs, cube_filename)
    
    # Show events
    show_events(events, mod, e_min, e_max, skydir, duration)
