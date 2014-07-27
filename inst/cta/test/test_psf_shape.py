#! /usr/bin/env python
# ==========================================================================
# This script simulates the PSF distribution.
#
# If matplotlib is installed, the PSF will be displayed on the screen.
#
# --------------------------------------------------------------------------
#
# Copyright (C) 2013-2014 Juergen Knoedlseder
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
import gammalib
import math
# Check if matplotlib is available
try:
    import matplotlib.pyplot as plt
    has_matplotlib = True
except ImportError:
    print("Matplotlib is not (correctly) installed on your system.")
    has_matplotlib = False


# ================ #
# Show Psf profile #
# ================ #
def show_psf(response, energy, r_max=0.8, rbins=1000):
    """
    Show Psf profile.
    """

    # Continue only if matplotlib is available
    if has_matplotlib:

        # Create figure
        plt.figure(1)
        plt.title("Point spread function ("+str(energy)+" TeV)")

        # Set log10(energy)
        logE = math.log10(energy)

        # Create delta axis
        delta = []
        dr    = r_max/rbins
        for i in range(rbins):
            #delta.append((i+0.5)*dr)
            delta.append((i)*dr)

        # Set Crab direction
        dir = gammalib.GSkyDir()
        dir.radec_deg(83.6331, 22.0145)
        eng = gammalib.GEnergy(energy, "TeV")

        # Get expected PSF
        sum       = 0.0
        psf       = []
        for i in range(rbins):
            r     = delta[i]*math.pi/180.0
            value = rsp.psf()(dir, r, eng) * \
                    2.0*math.pi * math.sin(r) * dr * math.pi/180.0
            sum += value
            psf.append(value)
        print(sum)

        # Plot PSF
        plt.plot(delta, psf, 'b-')
        plt.xlim([0.0,0.01])
        #plt.ylim([0.0,1.0e-18])
        
        # Set axes
        plt.xlabel("Offset angle (degrees)")
        plt.ylabel("Psf integral")

        # Notify
        print("PLEASE CLOSE WINDOW TO CONTINUE ...")

        # Show plot
        plt.show()

    # Return
    return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    Show Psf profile
    """
    # Dump header
    print("")
    print("********************")
    print("* Show Psf profile *")
    print("********************")

    # Load response
    #rsp = GCTAResponseIrf("irf_file.fits",
    #                      GCaldb("../caldb/data/cta/e/bcf/IFAE20120510_50h_King"))
    #rsp = GCTAResponseIrf("irf_file.fits",
    #                      GCaldb("../caldb/data/cta/e/bcf/IFAE20120510_50h"))
    #rsp = GCTAResponseIrf("cta_dummy_irf",
    #                      GCaldb("../caldb"))
    exposure = gammalib.GCTAExposure("data/expcube.fits")
    #psf      = gammalib.GCTAMeanPsf("data/psfcube.fits")
    psf      = gammalib.GCTAMeanPsf("psfcube.fits")
    rsp      = gammalib.GCTAResponseCube(exposure, psf)

    # Show PSF
    show_psf(rsp, 1.0)
    #show_psf(rsp, 1.0, r_max=0.2)
