#! /usr/bin/env python
# ==========================================================================
# This script tests the Npred computation of the CTA response.
#
# Requires:
# - matplotlib (optional)
#
# Copyright (C) 2011-2012 Jurgen Knodlseder
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
import math
import os
import glob
import sys

# Optionally import matplotlib if available
try:
    import matplotlib.pyplot as plt
    has_matplotlib = True
except ImportError:
    print("WARNING: Matplotlib not available for plotting.")
    has_matplotlib = False


# ======================================== #
# Compute NPsf as function of PSF position #
# ======================================== #
def npsf(caldb, irf, radius=1.0, energy=0.5):
    """
    Compute NPsf as function of PSF position.
    """
    # Setup CTA response. If no calibration database is given we directly
    # load the PSF from the file specified using the irf parameter.
    if caldb == "":
        rsp = GCTAResponse()
        rsp.load_psf(irf)
    else:
        rsp = GCTAResponse(irf, caldb)

    # Setup ROI
    instDir = GCTAInstDir()
    instDir.radec_deg(0.0, 0.0)
    roi = GCTARoi()
    roi.centre(instDir)
    roi.radius(radius)

    # Setup source energy
    srcEng = GEnergy()
    srcEng.TeV(energy)

    # Setup dummy time and pointing
    srcTime = GTime()
    pnt = GCTAPointing()

    # Compute PSF as function of offset
    x = []
    y = []
    rmax = int(radius * 150.0)
    for i in range(rmax):

        # Compute offset angle
        offset = i * radius / 100.0

        # Set PSF position
        srcDir = GSkyDir()
        srcDir.radec_deg(0.0, offset)

        # Compute NPsf
        npsf = rsp.npsf(srcDir, srcEng.log10TeV(), srcTime, pnt, roi)

        # Store results
        x.append(offset)
        y.append(npsf)
        # print offset, npsf

    # Optionally plot results
    if has_matplotlib:

        # Create figure
        plt.figure(1)
        plt.title("NPsf")

        # Plot data
        plt.plot(x, y, 'r-')
        plt.plot([radius, radius], [0.0, 1.0], 'b-')

        # Set axes
        plt.xlabel("Offset angle (deg)")
        plt.ylabel("Npsf")

        # Show plot
        plt.show()

    # Return
    return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    Test NPsf computation
    """
    # Set response information
    caldb = "$GAMMALIB/share/caldb/cta"
    irf   = "kb_E_50h_v3"

    # Compute NPsf
    npsf(caldb, irf, radius=1.0)
    # npsf("", file, radius=1.0, energy=1.0)
