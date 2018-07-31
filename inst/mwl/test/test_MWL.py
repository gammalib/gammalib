# ==========================================================================
# This module performs unit tests for the GammaLib MWL module.
#
# Copyright (C) 2012-2018 Juergen Knoedlseder
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
import gammalib
import test_support


# ================================== #
# Test class for GammaLib MWL module #
# ================================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for GammaLib MWL module
    """
    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Return
        return

    # Test class pickeling
    def _test_pickeling(self):
        """
        Test class pickeling
        """
        # Perform pickeling tests of empty classes
        test_support.pickeling(self, gammalib.GMWLDatum())
        test_support.pickeling(self, gammalib.GMWLObservation())
        test_support.pickeling(self, gammalib.GMWLSpectrum())

        # Setup test
        eng   = gammalib.GEnergy(1.0, 'TeV')
        datum = gammalib.GMWLDatum(eng, eng, 1.0, 0.1)
        obs   = gammalib.GMWLObservation(os.environ['TEST_MWL_DATA']+'/crab_mwl.fits')
        spec  = gammalib.GMWLSpectrum(os.environ['TEST_MWL_DATA']+'/crab_mwl.fits')

        # Perform pickeling tests of filled classes
        test_support.pickeling(self, gammalib.GMWLDatum(datum))
        test_support.pickeling(self, gammalib.GMWLObservation(obs))
        test_support.pickeling(self, gammalib.GMWLSpectrum(spec))

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('MWL')

        # Append tests
        self.append(self._test_pickeling, 'Test MWL class pickeling')

        # Return
        return
