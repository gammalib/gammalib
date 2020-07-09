#! /usr/bin/env python
# ==========================================================================
# Compare radial model to MC
# --------------------------------------------------------------------------
#
# Copyright (C) 2020 Juergen Knoedlseder
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


# ========================== #
# Compare radial model to MC #
# ========================== #
def compare(radial, rad_max=1.0, nrad=100, nmc=1000000):
    """
    Compare radial model to MC.
    """
    # Get dummy energy and time
    energy = gammalib.GEnergy()
    time   = gammalib.GTime()
    
    # Generate theta vector
    theta_bin = rad_max/nrad
    theta     = [(i+0.5) * theta_bin for i in range(nrad)]

    # Generate model vector
    model = [radial.eval(t*gammalib.deg2rad, energy, time) for t in theta]

    # Generate MC vector
    ran = gammalib.GRan()
    mc  = [0.0 for t in theta]
    for i in range(nmc):
        dir    = radial.mc(energy, time, ran)
        offset = radial.dir().dist_deg(dir)
        index  = int(offset / rad_max * nrad)
        if index >= 0 and index < nrad:
            mc[index] += 1.0

    # Divide MC by area of ring
    for i, t in enumerate(theta):
        tmin  = t-0.5*theta_bin
        tmax  = t+0.5*theta_bin
        area  = gammalib.pi * (tmax*tmax - tmin*tmin)
        mc[i] = mc[i] / area

    # Normalise model
    n_model = sum(model)
    n_mc    = sum(mc)
    if n_model > 0.0:
        norm = n_mc/n_model
        model = [m*norm for m in model]

    # Create figure
    plt.figure(1)
    plt.title('Model-MC comparison')

    # Plot data
    plt.plot(theta, model, 'r-')
    plt.plot(theta, mc,    'b-')

    # Set axes
    plt.xlabel('Radius (deg)')
    plt.ylabel('Number of photons')

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

    # Define model centre
    centre = gammalib.GSkyDir()
    centre.radec_deg(83.6331, 22.0145)

    # Define radial model
    radial = gammalib.GModelSpatialRadialRing(centre, 0.45, 0.15)

    # Compare model to MC
    compare(radial)
