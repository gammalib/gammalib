# ==========================================================================
# This module performs unit tests for the GammaLib CTA module.
#
# Copyright (C) 2012-2014 Juergen Knoedlseder
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


# ================================== #
# Test class for GammaLib CTA module #
# ================================== #
class Test(GPythonTestSuite):
    """
    Test class for GammaLib CTA module.
    """
    # Constructor
    def __init__(self):
        """
        Constructor.
        """
        # Call base class constructor
        GPythonTestSuite.__init__(self)

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name("CTA")

        # Append tests
        self.append(self.test_aeff, "Test CTA effective area classes")
        self.append(self.test_psf, "Test CTA PSF classes")
        self.append(self.test_edisp, "Test CTA energy dispersion classes")
        self.append(self.test_onoff, "Test CTA ON/OFF analysis")

        # Return
        return

    # Test effective area response
    def test_aeff(self):
        """
        Test GCTAAeff classes.
        """
        # Load 2D Aeff
        self.test_try("Test GCTAAeff2D file constructor")
        try:
            db = GCaldb("../inst/cta/caldb")
            db.open("cta", "e")
            filename = db.filename("","","EFF_AREA","","","NAME(IFAE20120510_50h)")
            aeff     = GCTAAeff2D(filename)
            self.test_try_success()
        except:
            self.test_try_failure("Unable to allocate GCTAAeff2D from file.")

        # Test Aeff values
        self.test_value(aeff(0.0, 0.0), 5535774176.75, 0.1)
        self.test_value(aeff(1.0, 0.0), 20732069462.7, 0.1)
        self.test_value(aeff(0.0, 0.01745), 5682897797.76, 0.1)
        self.test_value(aeff(1.0, 0.01745), 18446656815.1, 0.1)

        # Load performance file
        self.test_try("Test GCTAAeffPerfTable file constructor")
        try:
            filename = "../inst/cta/caldb/cta_dummy_irf.dat"
            aeff     = GCTAAeffPerfTable(filename)
            self.test_try_success()
        except:
            self.test_try_failure("Unable to allocate GCTAAeffPerfTable from file.")

        # Test Aeff values
        self.test_value(aeff(0.0, 0.0), 2738898000.0, 0.1)
        self.test_value(aeff(1.0, 0.0), 16742420500.0, 0.1)
        self.test_value(aeff(0.0, 0.01745), 2590995083.29, 0.1)
        self.test_value(aeff(1.0, 0.01745), 15838314971.2, 0.1)

        # Load 1DC ARF response file
        self.test_try("Test GCTAAeffArf file constructor")
        try:
            filename = "../inst/cta/test/caldb/dc1/arf.fits"
            aeff     = GCTAAeffArf(filename)
            self.test_try_success()
        except:
            self.test_try_failure("Unable to allocate GCTAAeffArf from file.")

        # Test Aeff values
        self.test_value(aeff(0.0, 0.0), 1607246236.98, 0.1)
        self.test_value(aeff(1.0, 0.0), 4582282342.98, 0.1)
        self.test_value(aeff(0.0, 0.01745), 1607246236.98, 0.1)
        self.test_value(aeff(1.0, 0.01745), 4582282342.98, 0.1)

        # Return
        return

    # Test point spread function response
    def test_psf(self):
        """
        Test GCTAPsf classes.
        """
        # Load 2D PSF
        self.test_try("Test GCTAPsf2D file constructor")
        try:
            db = GCaldb("../inst/cta/caldb")
            db.open("cta", "e")
            filename = db.filename("","","RPSF","","","NAME(IFAE20120510_50h)")
            psf      = GCTAPsf2D(filename)
            self.test_try_success()
        except:
            self.test_try_failure("Unable to allocate GCTAPsf2D from file.")

        # Test PSF values
        self.test_value(psf(0.0, 0.0, 0.0), 71870.337730, 1.0e-6)
        self.test_value(psf(0.001, 0.0, 0.0), 57344.463816, 1.0e-6)
        self.test_value(psf(0.0, 1.0, 0.0), 270344.329947, 1.0e-6)
        self.test_value(psf(0.001, 1.0, 0.0), 115628.756239, 1.0e-6)
        self.test_value(psf(0.0, 1.0, 0.01745), 190115.769257, 1.0e-6)
        self.test_value(psf(0.001, 1.0, 0.01745), 104623.364268, 1.0e-6)

        # Load King profile PSF
        self.test_try("Test GCTAPsfKing file constructor")
        try:
            db = GCaldb("../inst/cta/caldb")
            db.open("cta", "e")
            filename = db.filename("","","RPSF","","","NAME(IFAE20120510_50h_King)")
            psf      = GCTAPsfKing(filename)
            self.test_try_success()
        except:
            self.test_try_failure("Unable to allocate GCTAPsfKing from file.")

        # Test PSF values
        self.test_value(psf(0.0, 0.0, 0.0), 213666.253408, 1.0e-6)
        self.test_value(psf(0.001, 0.0, 0.0), 90939.5585827, 1.0e-6)
        self.test_value(psf(0.0, 1.0, 0.0), 1127688.29309, 1.0e-5)
        self.test_value(psf(0.001, 1.0, 0.0), 54916.6796886, 1.0e-6)
        self.test_value(psf(0.0, 1.0, 0.01745), 660973.856068, 1.0e-6)
        self.test_value(psf(0.001, 1.0, 0.01745), 80272.5530008, 1.0e-6)

        # Load performance file
        self.test_try("Test GCTAPsfPerfTable file constructor")
        try:
            filename = "../inst/cta/caldb/cta_dummy_irf.dat"
            psf      = GCTAPsfPerfTable(filename)
            self.test_try_success()
        except:
            self.test_try_failure("Unable to allocate GCTAPsfPerfTable from file.")

        # Test PSF values
        self.test_value(psf(0.0, 0.0, 0.0), 537855.359317, 1.0e-6)
        self.test_value(psf(0.001, 0.0, 0.0), 99272.364544, 1.0e-6)
        self.test_value(psf(0.0, 1.0, 0.0), 1292609.56448, 1.0e-6)
        self.test_value(psf(0.001, 1.0, 0.0), 22277.2429186, 1.0e-6)
        self.test_value(psf(0.0, 1.0, 0.01745), 1292609.56448, 1.0e-6)
        self.test_value(psf(0.001, 1.0, 0.01745), 22277.2429186, 1.0e-6)

        # Load 1DC PSF file
        self.test_try("Test GCTAPsfVector file constructor")
        try:
            filename = "../inst/cta/test/caldb/dc1/psf_magic.fits"
            psf      = GCTAPsfVector(filename)
            self.test_try_success()
        except:
            self.test_try_failure("Unable to allocate GCTAPsfVector from file.")

        # Print
        self.test_value(psf(0.0, -1.0, 0.0), 42264.1147434, 1.0e-6)
        self.test_value(psf(0.001, -1.0, 0.0), 37009.0228006, 1.0e-6)
        self.test_value(psf(0.0, 0.0, 0.0), 208989.943127, 1.0e-6)
        self.test_value(psf(0.001, 0.0, 0.0), 108388.810748, 1.0e-6)
        self.test_value(psf(0.0, 0.0, 0.01745), 208989.943127, 1.0e-6)
        self.test_value(psf(0.001, 0.0, 0.01745), 108388.810748, 1.0e-6)

        # Return
        return

    # Test point spread function response
    def test_edisp(self):
        """
        Test GCTAEdisp classes.
        """
        # Load performance file
        self.test_try("Test GCTAEdispPerfTable file constructor")
        try:
            edisp = GCTAEdispPerfTable("../inst/cta/test/caldb/cta_dummy_irf.dat")
            self.test_try_success()
        except:
            self.test_try_failure("Unable to allocate GCTAEdispPerfTable from file.")

        # Test energy dispersion values
        self.test_value(edisp(0.0, 0.0), 9.99019627861, 1.0e-6)
        self.test_value(edisp(0.001, 0.0), 9.9870644077, 1.0e-6)
        self.test_value(edisp(0.01, 0.0), 9.68182, 1.0e-6)
        self.test_value(edisp(0.1, 0.0), 0.434382, 1.0e-6)
        self.test_value(edisp(1.0, 1.0), 18.064868197, 1.0e-6)
        self.test_value(edisp(1.001, 1.0, 0.0), 18.0463571212, 1.0e-6)

        # Load response
        self.test_try("Test GCTAResponse file constructor")
        try:
            rsp = GCTAResponse("cta_dummy_irf", GCaldb("../inst/cta/test/caldb"))
            self.test_try_success()
        except:
            self.test_try_failure("Unable to allocate GCTAResponse from file.")

        # Test nedisp computations
        dir  = GSkyDir()
        pnt  = GCTAPointing()
        time = GTime()
        self.test_value(rsp.nedisp(dir, GEnergy(3.7, "TeV"), time, pnt, 
                                   GEbounds(GEnergy(0.1, "TeV"),
                                            GEnergy(10.0, "TeV"))),
                        1.0, 1.0e-6)
        self.test_value(rsp.nedisp(dir, GEnergy(3.7, "TeV"), time, pnt, 
                                   GEbounds(GEnergy(2.72345,  "TeV"),
                                            GEnergy(5.026615, "TeV"))),
                        1.0, 1.0e-6)
        self.test_value(rsp.nedisp(dir, GEnergy(3.7, "TeV"), time, pnt, 
                                   GEbounds(GEnergy(3.7, "TeV"),
                                            GEnergy(10.0, "TeV"))),
                        0.5, 1.0e-6)

    # Test ON/OFF analysis
    def test_onoff(self):
        """
        Test ON/OFF analysis.
        """
        # Create ON region
        crab = GSkyDir()
        crab.radec_deg(83.6331, 22.0145)
        on = GSkyRegions()
        on.append(GSkyRegionCircle(crab, 1.0))

        # Create OFF region
        off = on
        
        # Set energy binning
        etrue = GEbounds(10, GEnergy(0.1, "TeV"), GEnergy(10.0, "TeV"))
        ereco = GEbounds(10, GEnergy(0.1, "TeV"), GEnergy(10.0, "TeV"))
        
        # Create ON/OFF observations by filling all events found in
        # the observation container and computing the response
        obs    = GObservations("../inst/cta/test/data/obs_unbinned.xml")
        onoffs = GCTAOnOffObservations()
        self.test_try("Test ON/OFF observation creation")
        try:
            for run in obs:
                onoff = GCTAOnOffObservation(ereco, on, off)
                onoff.fill(run)
                onoff.compute_response(run, etrue)
                onoffs.append(onoff)
            self.test_try_success()
        except:
            self.test_try_failure("Unable to create ON/OFF observations.")

        # Save PHA, ARF and RMFs
        self.test_try("Test ON/OFF observation PHA, ARF and RMF saving")
        try:
            for run in onoffs:
                run.on_spec().save("onoff_pha_on.fits", True)
                run.off_spec().save("onoff_pha_off.fits", True)
                run.arf().save("onoff_arf.fits", True)
                run.rmf().save("onoff_rmf.fits", True)
            self.test_try_success()
        except:
            self.test_try_failure("Unable to save ON/OFF PHA, ARF and RMFs.")

        # Save ON/OFF observations
        self.test_try("Test ON/OFF observation saving")
        try:
            onoffs.save("onoff.xml")
            self.test_try_success()
        except:
            self.test_try_failure("Unable to save ON/OFF observations.")
        
        # Return
        return
