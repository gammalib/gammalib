from gammalib import *
from ctools import *
import math
import os
import matplotlib.pyplot as plt
import numpy


# def LiMa(n_observed, mu_background):
#     term_a = math.sign(n_observed - mu_background) * sqrt(2)
#     term_b = math.sqrt(n_observed * log(n_observed / mu_background) - n_observed + mu_background)
#     return term_a * term_b
    
def searchctr(deltasqs,psfvals):
    ctr = 0
    sum_total = sum(psfvals)
    for i in range(len(psfvals)):
        sum1 = sum(psfvals[0:i])
        if sum1 > sum_total*0.68:
            ctr = math.sqrt(deltasqs[i])
            break
    return ctr

# ============= #
# Plot spectrum #
# ============= #
def plot_spectrum(observations, emin, emax, specbins=10):
    """
    Plot the observed and the fitted counts spectra.

    Parameters:
     observations - Observations container
    Keywords:
     emin         - Lower event energy (TeV)
     emax         - Largest event energy (TeV)
     specbins     - Number of spectral bins
    """
    # Set legend fontsize
    params = {'legend.fontsize': 10}
    plt.rcParams.update(params)

    # Set plot styles
    styles = ['b-', 'g-', 'y-', 'n-']

    # Set plot title
    plt.title("Spectrum")

    # Create logarithmic energy axis in TeV
    ebounds = GEbounds()
    e_min = GEnergy()
    e_max = GEnergy()
    e_min.TeV(emin)
    e_max.TeV(emax)
    ebounds.set_log(specbins, e_min, e_max)
    energy = [ebounds.elogmean(i).TeV() for i in range(ebounds.size())]

    # Initialise spectra
    nmodels = observations.models().size()
    counts = [0.0 for i in range(ebounds.size())]
    model = [0.0 for i in range(ebounds.size() * nmodels)]
    sum_model = [0.0 for i in range(ebounds.size())]

    # Loop over observations
    for obs in observations:

        # Make sure that we have a CTA observation
        obs = GCTAObservation(obs)

        # Get event list
        list = obs.events()

        # Get livetime
        time = GTime()  # Dummy time needed for npred
        livetime = obs.livetime()

        # Set counters
        sum_cts = 0.0
        sum_mod = 0.0

        # Create spectrum
        for atom in list:
            index = ebounds.index(atom.energy())
            counts[index] = counts[index] + 1.0
            sum_cts += 1.0

        # Create error bars
        error = [math.sqrt(c) for c in counts]

        # Extract models
        for k, m in enumerate(observations.models()):
            ioff = k * ebounds.size()
            for i in range(ebounds.size()):
                eng = ebounds.elogmean(i)
                ebinsize = ebounds.emax(i) - ebounds.emin(i)
                expect         = m.npred(eng, time, obs) * \
                    ebinsize.MeV() * livetime
                model[i + ioff] += expect
                sum_mod += expect

        # Print results
        print("-", sum_cts, sum_mod)

    # Create sum spectrum
    for k, m in enumerate(observations.models()):
        ioff = k * ebounds.size()
        for i in range(ebounds.size()):
            sum_model[i] += model[i + ioff]

    # Plot spectrum
    plt.loglog(energy, counts, 'ro', label='data')
    plt.errorbar(energy, counts, error, fmt=None, ecolor='r')
    for k, m in enumerate(observations.models()):
        ioff = k * ebounds.size()
        plt.loglog(energy, model[ioff:ioff + ebounds.size()],
                   styles[k], label=m.name())
    plt.loglog(energy, sum_model, 'r-', label='total')

    # Put labels
    plt.xlabel("Energy (TeV)")
    plt.ylabel("Counts")
    plt.legend(loc="lower left")

    # Return
    return


# =========== #
# Plot offset #
# =========== #
def plot_offset(observations, emin, emax, ra, dec, rad, offsetbins=30):
    """
    Plot observed and fitted offset angles.

    Parameters:
     observations - Observations container
     emin         - Lower event energy (TeV)
     emax         - Largest event energy (TeV)
     ra           - Right Ascension of ROI centre (deg)
     dec          - Declination of ROI centre (deg)
     rad          - Radius of ROI (deg)
 offsetbins   - Number of bins for offset histogram
    """
    # Set legend fontsize
    params = {'legend.fontsize': 10}
    plt.rcParams.update(params)

    # Set plot styles
    styles = ['b-', 'g-', 'y-', 'n-']

    # Set plot title
    plt.title("Offset")

    # Set binsize
    binsz = rad / offsetbins
    if binsz < 0.01:
        binsz = 0.01

    # Set parameters
    nxpix = int(2.0 * rad / binsz)
    nypix = int(2.0 * rad / binsz)
    enumbins = 10
    coordsys = "CEL"
    proj = "CAR"

    # Create offset angle distribution
    nmodels = observations.models().size()
    doffset = rad * rad / offsetbins
    offset = [(i + 0.5) * doffset for i in range(offsetbins)]
    counts = [0.0 for i in range(offsetbins)]
    model = [0.0 for i in range(offsetbins * nmodels)]
    sum_model = [0.0 for i in range(offsetbins)]

    # Loop over observations
    for iobs, obs in enumerate(observations):

        # Get CTA event cube
        cube = GCTAEventCube(obs.events())

        # Get pointing direction
        pnt = obs.pointing().dir()

        # Get ROI centre and radius
        roi_dir = GSkyDir()
        roi_dir.radec_deg(ra, dec)
        roi_rad = rad

        # Create offset histogram with respect to ROI centre
        sum_cts = 0.0
        sum_mod = 0.0
        for bin in cube:
            off = bin.dir().dist_deg(roi_dir)
            if off < rad:
                off2 = off * off
                index = int(off2 / doffset)
                if index < offsetbins:
                    counts[index] += bin.counts()
                    sum_cts += bin.counts()
                    for k, m in enumerate(observations.models()):
                        model[index + k *
                              offsetbins] += m.eval(bin, obs) * bin.size()
                        sum_mod += m.eval(bin, obs) * bin.size()
        print("-", iobs, sum_cts, sum_mod)

    # Create sum histogram
    for k, m in enumerate(observations.models()):
        ioff = k * offsetbins
        for i in range(offsetbins):
            sum_model[i] += model[i + ioff]

    # Create error bars
    error = [math.sqrt(c) for c in counts]
    # Plot
    plt.plot(offset, counts, 'ro', label='data')
    plt.errorbar(offset, counts, error, fmt=None, ecolor='r')
    for k, m in enumerate(observations.models()):
        plt.plot(offset, model[k * offsetbins:k * offsetbins + offsetbins],
                 styles[k], label=m.name())
    plt.plot(offset, sum_model, 'r-', label='total')

    # Put labels
    plt.xlabel("Offset squared (deg^2)")
    plt.ylabel("Counts")
    plt.legend(loc="upper right")

    # Return
    return


# ======================== #
# Plot residual counts map #
# ======================== #
def plot_residuals(observations, emin, emax, ra, dec, rad):
    """
    Plot residual counts map.

    Parameters:
     observations - Observations container
     emin         - Lower event energy (TeV)
     emax         - Largest event energy (TeV)
     ra           - Right Ascension of ROI centre (deg)
     dec          - Declination of ROI centre (deg)
     rad          - Radius of ROI (deg)
    """
 
    # Set binsize
    binsz = rad / 30.0
    if binsz < 0.01:
        binsz = 0.01

    # Set ctbin parameters
    nxpix = int(2.0 * rad / binsz)
    nypix = int(2.0 * rad / binsz)
    enumbins = 10
    coordsys = "CEL"
    proj = "CAR"


    # Initialise maps
    counts = numpy.empty((nxpix, nypix), dtype=numpy.float64) * 0.0
    model = numpy.empty((nxpix, nypix), dtype=numpy.float64) * 0.0
    residual = numpy.empty((nxpix, nypix), dtype=numpy.float64) * 0.0
    significance = numpy.empty((nxpix, nypix), dtype=numpy.float64) * 0.0

    # Set map boundaries for display
    rabin = binsz / math.cos(dec * math.pi / 180.0)
    decbin = binsz
    dra = rabin * nxpix
    ddec = decbin * nypix
    ra_min = ra - dra / 2.0
    ra_max = ra + dra / 2.0
    dec_min = dec - ddec / 2.0
    dec_max = dec + ddec / 2.0

    # Loop over observations
    for iobs, obs in enumerate(observations):

        # Make sure that we have a CTA observation
        obs = GCTAObservation(obs)

        # Get counts map
        events = GCTAEventCube(obs.events())
        cntmap = events.map().copy()

        # Get ROI centre and radius
        roi_dir = GSkyDir()
        roi_dir.radec_deg(ra, dec)
        roi_rad = rad

        # Make model map
        for event in events:
            sep = event.dir().dir().dist_deg(roi_dir)
            if sep <= roi_rad:
                mod = observations.models().eval(event, obs) * event.size()
                event.counts(mod)
            else:
                event.counts(0.0)

        # Fill the maps
        dir = GSkyDir()
        for ix in range(nxpix):
            ra_val = ra_max - ix * rabin
            for iy in range(nypix):
                dec_val = dec_max - iy * decbin
                dir.radec_deg(ra_val, dec_val)
                ipix = cntmap.dir2pix(dir)
                index = cntmap.dir2inx(dir)
                if ipix.x() >= 0 and ipix.x() < cntmap.nx() and \
                   ipix.y() >= 0 and ipix.y() < cntmap.ny() and index < cntmap.npix():
                    for iebin in range(enumbins):
                        counts[iy, ix] += cntmap[index, iebin]
                        model[iy, ix] += events.map()[index, iebin]
        print("-", iobs, counts.sum(), model.sum())

    # Create residual map
    #residual = counts - model
    print("Residual: ", counts.sum() - model.sum())
    for ix in range(nxpix):
        for iy in range(nypix):
            if model[iy, ix] > 0.0:
                residual[iy, ix] = (counts[iy, ix] - model[iy, ix]) / \
                    math.sqrt(model[iy, ix])

    # Create significance map
    print("Creating significance map using LiMa formula")
  #  for ix in range(nxpix):
   #     for iy in range(nypix):
    #        significance[iy, ix] = LiMa(counts[iy, ix], model[iy,ix])

    # Show counts map
    plt.subplot(221)
    plt.title("Counts map")
    im = plt.imshow(counts, interpolation='nearest',
                    extent=[ra_max, ra_min, dec_min, dec_max])
    plt.colorbar(im, orientation='horizontal', shrink=0.8)

    # Put labels
    plt.xlabel("Right Ascension (deg)")
    plt.ylabel("Declination (deg)")

    # Show model map
    plt.subplot(222)
    plt.title("Model map")
    im = plt.imshow(model, interpolation='nearest',
                    extent=[ra_max, ra_min, dec_min, dec_max])
    plt.colorbar(im, orientation='horizontal', shrink=0.8)

    # Put labels
    plt.xlabel("Right Ascension (deg)")
    plt.ylabel("Declination (deg)")


    # Show residual map
    plt.subplot(224)
    plt.title("Residual map (rms)")
    im = plt.imshow(residual, interpolation='nearest',
                    extent=[ra_max, ra_min, dec_min, dec_max])
    plt.colorbar(im, orientation='horizontal', shrink=0.8)

    # Put labels
    plt.xlabel("Right Ascension (deg)")
    plt.ylabel("Declination (deg)")

    # Return
    return


