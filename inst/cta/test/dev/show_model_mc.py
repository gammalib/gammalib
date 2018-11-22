#! /usr/bin/env python
# ==========================================================================
# Show simulation - model for a given model
#
# Copyright (C) 2018 Juergen Knoedlseder
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
import ctools
import cscripts
from cscripts import obsutils
import matplotlib.pyplot as plt


# ====================== #
# Plot residual spectrum #
# ====================== #
def plot_resspec(frame, energies, resspec, errspec, emin=0.1, emax=100.0):
    """
    Plot residual spectrum

    Parameters
    ----------
    frame : pyplot
        Frame for spectrum
    """
    # Set axes scales and limit
    frame.set_xscale('log')
    frame.set_xlim([emin, emax])

    # Counts and model
    frame.errorbar(energies, resspec, yerr=errspec,
                   fmt='ko', capsize=0, linewidth=2, zorder=2, label='Data')
    frame.plot(frame.get_xlim(),[0,0], 'r-')
    frame.set_xlabel('Energy (TeV)')
    frame.set_ylabel('Simulation - Model')

    # Return
    return


# ================= #
# Plot residual map #
# ================= #
def plot_resmap(plt1, plt2, map, error, steps=1):
    """
    Plot pull histogram

    Parameters
    ----------
    plt1 : pyplot
        Frame for first map
    plt2 : pyplot
        Frame for second map
    map : `~gammalib.GSkyMap()`
        Sky map
    error : `~gammalib.GSkyMap()`
        Sky map
    steps : int
        Number of steps for rebinning
    """
    # Create longitude array from skymap
    x_lon = []
    y_lon = []
    e_lon = []
    for ix in range(0, map.nx(), steps):
        value = 0.0
        err   = 0.0
        for istep in range(steps):
            if istep+ix < map.nx():
                for iy in range(map.ny()):
                    index  = istep+ix+iy*map.nx()
                    value += map[index]
                    err   += error[index]*error[index]
        x_lon.append(ix)
        y_lon.append(value)
        e_lon.append(math.sqrt(err))

    # Create latitude array from skymap
    x_lat = []
    y_lat = []
    e_lat = []
    for iy in range(0, map.ny(), steps):
        value = 0.0
        err   = 0.0
        for istep in range(steps):
            if istep+iy < map.ny():
                for ix in range(map.nx()):
                    index  = ix+(istep+iy)*map.nx()
                    value += map[index]
                    err   += error[index]*error[index]
        x_lat.append(iy)
        y_lat.append(value)
        e_lat.append(math.sqrt(err))

    # Plot longitude distribution
    plt1.errorbar(x_lon, y_lon, yerr=e_lon,
                  fmt='ko', capsize=0, linewidth=2, zorder=2)
    plt1.axhline(0, color='r', linestyle='--')
    plt1.set_xlabel('Right Ascension (pixels)')
    plt1.set_ylabel('Simulation - Model')
    plt1.grid(True)

    # Plot latitude distribution
    plt2.errorbar(x_lat, y_lat, yerr=e_lat,
                  fmt='ko', capsize=0, linewidth=2, zorder=2)
    plt2.axhline(0, color='r', linestyle='--')
    plt2.set_xlabel('Declination (pixels)')
    plt2.set_ylabel('Simulation - Model')
    plt2.grid(True)

    # Return
    return


# ======== #
# Plot map #
# ======== #
def plot_map(frame, map, smooth=0.0):
    """
    Plot Aitoff map

    Parameters
    ----------
    frame : pyplot
        Frame for map
    map : `~gammalib.GSkyMap()`
        Sky map
    smooth : float, optional
        Map smoothing parameter (degrees)
    """
    # Optionally smooth map
    if smooth > 0.0:
        _map = map.copy()
        _map.smooth('DISK', smooth)
    else:
        _map = map
    
    # Create array from skymap
    array = []
    v_max = 0.0
    for iy in range(_map.ny()):
        row = []
        for ix in range(_map.nx()):
            index = ix+iy*_map.nx()
            value = _map[index]
            row.append(value)
        array.append(row)

    # Show Aitoff projection
    c = frame.imshow(array, aspect=1.0)
    cbar = plt.colorbar(c, orientation='vertical', shrink=0.7)
    frame.grid(True)
    frame.set_xlabel('Right Ascension (pixels)')
    frame.set_ylabel('Declination (pixels)')

    # Return
    return


# =========================== #
# Determine spatial residuals #
# =========================== #
def residual_spatial(countscube, modelcube):
    """
    Determine spatial residuals

    Parameters
    ----------
    countscube : `~gammalib.GCTAEventCube`
        Simulated counts cube
    modelcube : `~gammalib.GCTAEventCube`
        Model cube

    Returns
    -------
    map, error : `~gammalib.GSkyMap()`
        Residual sky map and error
    """
    # Get model cube maps
    map1 = countscube.counts()
    map2 = modelcube.counts()

    # Generate residual map
    map = map1 - map2
    map.stack_maps()

    # Generate residual error bars
    error = map2.copy()
    error.stack_maps()
    error = error.sqrt()

    # Return residual map and error map
    return map, error


# ============================ #
# Determine spectral residuals #
# ============================ #
def residual_spectral(countscube, modelcube):
    """
    Determine spectral residuals

    Parameters
    ----------
    countscube : `~gammalib.GCTAEventCube`
        Simulated counts cube
    modelcube : `~gammalib.GCTAEventCube`
        Model cube

    Returns
    -------
    map, error : list
        Residual vector
    """
    # Get model cube maps
    map1 = countscube.counts()
    map2 = modelcube.counts()

    # Compute spectra
    spec1    = []
    spec2    = []
    residual = []
    error    = []
    energies = []
    for ieng in range(map1.nmaps()):
        sum1 = 0.0
        sum2 = 0.0
        for ipix in range(map1.npix()):
            sum1 += map1[ipix, ieng]
            sum2 += map2[ipix, ieng]
        spec1.append(sum1)
        spec2.append(sum2)
        residual.append(sum1 - sum2)
        error.append(math.sqrt(sum2))
        energies.append(modelcube.energy(ieng).TeV())

    # Return residual and error
    return energies, residual, error


# ==================== #
# Simulate counts cube #
# ==================== #
def simulate_counts_cube(obs, ra=83.63, dec=22.01, npix=50, binsz=0.02,
                         emin=0.1, emax=100.0, ebins=20, seed=1, edisp=False):
    """
    Simulate counts cube

    Parameters
    ----------
    obs : `~gammalib.GObservations`
        Observation container
    ebins : int, optional
        Number of energy bins for binned or stacked analysis
    edisp : boolean, optional
        Use energy dispersion
    """
    # Simulate events
    ctobssim = ctools.ctobssim(obs)
    ctobssim['seed']  = seed
    ctobssim['edisp'] = edisp
    ctobssim.run()
    
    # Create counts cube
    ctbin = ctools.ctbin(ctobssim.obs())
    ctbin['stack']    = True
    ctbin['ebinalg']  = 'LOG'
    ctbin['emin']     = emin
    ctbin['emax']     = emax
    ctbin['enumbins'] = ebins
    ctbin['coordsys'] = 'CEL'
    ctbin['proj']     = 'CAR'
    ctbin['xref']     = ra
    ctbin['yref']     = dec
    ctbin['nxpix']    = npix
    ctbin['nypix']    = npix
    ctbin['binsz']    = binsz
    ctbin.run()

    # Return counts cube
    return ctbin.cube().copy()


# ================= #
# Create model cube #
# ================= #
def create_model_cube(obs, countscube, edisp=False):
    """
    Create model cube

    Parameters
    ----------
    obs : `~gammalib.GObservations`
        Observation container
    countscube : `~gammalib.GCTAEventCube`
        Counts cube for binning
    edisp : boolean, optional
        Use energy dispersion
    """
    # Create model cube
    ctmodel = ctools.ctmodel(obs)
    ctmodel.cube(countscube)
    ctmodel['edisp'] = edisp
    ctmodel.run()

    # Return model cube
    return ctmodel.cube().copy()


# ================= #
# Show model and MC #
# ================= #
def show_model_mc(modelname, duration=180000.0, caldb='prod2', irf='South_50h',
                  ebins=20, edisp=False):
    """
    """
    # Set pointing direction
    pntdir = gammalib.GSkyDir()
    pntdir.radec_deg(83.63, 22.01)
    
    # Set single CTA observation
    run = obsutils.set_obs(pntdir, duration=duration, caldb=caldb, irf=irf)

    # Set observation container
    obs = gammalib.GObservations()
    obs.append(run)

    # Load and append model
    models = gammalib.GModels(modelname)
    obs.models(models)

    # Create simulated counts cube
    countscube = simulate_counts_cube(obs, npix=400, binsz=0.002, ebins=ebins,
                                      edisp=edisp)
    print(countscube)

    # Create model cube
    modelcube = create_model_cube(obs, countscube, edisp=edisp)
    print(modelcube)

    # Create residual map and spectrum
    resmap, errmap             = residual_spatial(countscube, modelcube)
    energies, resspec, errspec = residual_spectral(countscube, modelcube)

    # Create figure
    fig = plt.figure(figsize=(14,5))

    # Create subplot for residual spectrum
    frame1 = fig.add_subplot(131)
    plot_resspec(frame1, energies, resspec, errspec, emin=0.05, emax=150.0)

    # Create subplots for residual profiles
    frame2 = fig.add_subplot(232)
    frame3 = fig.add_subplot(235)
    plot_resmap(frame2, frame3, resmap, errmap, steps=5)

    # Create subplot for residual map
    frame4 = fig.add_subplot(133)
    plot_map(frame4, resmap, smooth=0.02)

    # Show figure
    plt.show()

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Set model name
    modelname = 'crab.xml'

    # Show model and MC
    #show_model_mc(modelname, ebins=50, edisp=False)
    show_model_mc(modelname, ebins=50, edisp=True)
