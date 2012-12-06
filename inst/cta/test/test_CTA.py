# ==========================================================================
# This module performs unit tests for the GammaLib CTA module.
#
# Copyright (C) 2012 Juergen Knoedlseder
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
            aeff = GCTAAeff2D("../inst/cta/caldb/data/cta/e/bcf/000001/irf_test.fits")
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
            aeff = GCTAAeffPerfTable("../inst/cta/caldb/kb_E_50h_v3.dat")
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
            aeff = GCTAAeffArf("../inst/cta/test/caldb/dc1/arf.fits")
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
            psf = GCTAPsf2D("../inst/cta/caldb/data/cta/e/bcf/000001/irf_test.fits")
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

        # Load performance file
        self.test_try("Test GCTAPsfPerfTable file constructor")
        try:
            psf = GCTAPsfPerfTable("../inst/cta/caldb/kb_E_50h_v3.dat")
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
            psf = GCTAPsfVector("../inst/cta/test/caldb/dc1/psf_magic.fits")
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
