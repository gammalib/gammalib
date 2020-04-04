# ==========================================================================
# This module performs unit tests for the GammaLib INTEGRAL/SPI module.
#
# Copyright (C) 2020 by Juergen Knoedlseder
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
import test_support


# =========================================== #
# Test class for GammaLib INTEGRAL/SPI module #
# =========================================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for GammaLib INTEGRAL/SPI module
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

    # Dummy test function
    def _dummy_test(self):
        """
        Dummy test function
        """
        # Return
        return

    # Test class pickeling
    def _test_pickeling(self):
        """
        Test class pickeling
        """
        # Perform pickeling tests of empty classes
        test_support.pickeling(self, gammalib.GSPIEventBin())
        test_support.pickeling(self, gammalib.GSPIEventCube())
        test_support.pickeling(self, gammalib.GSPIInstDir())
        test_support.pickeling(self, gammalib.GSPIObservation())
        test_support.pickeling(self, gammalib.GSPIResponse())

        # Setup test (TODO: to be filled with meaningful values)
        bin  = gammalib.GSPIEventBin()
        cube = gammalib.GSPIEventCube()
        dir  = gammalib.GSPIInstDir()
        obs  = gammalib.GSPIObservation()
        rsp  = gammalib.GSPIResponse()

        # Perform pickeling tests of filled classes
        test_support.pickeling(self, gammalib.GSPIEventBin(bin))
        test_support.pickeling(self, gammalib.GSPIEventCube(cube))
        test_support.pickeling(self, gammalib.GSPIInstDir(dir))
        test_support.pickeling(self, gammalib.GSPIObservation(obs))
        test_support.pickeling(self, gammalib.GSPIResponse(rsp))

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('INTEGRAL/SPI')

        # Append tests
        self.append(self._dummy_test, 'INTEGRAL/SPI dummy test')
        self.append(self._test_pickeling, 'Test INTEGRAL/SPI class pickeling')

        # Return
        return
