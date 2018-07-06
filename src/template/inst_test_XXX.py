# ==========================================================================
# This module performs unit tests for the GammaLib [INSTRUMENT] module.
#
# Copyright (C) [YEAR] by [AUTHOR]
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


# ====================================== #
# Test class for GammaLib [INSTRUMENT] module #
# ====================================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for GammaLib [INSTRUMENT] module
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
        test_support._pickeling(self, gammalib.GXXXEventAtom())
        test_support._pickeling(self, gammalib.GXXXEventBin())
        test_support._pickeling(self, gammalib.GXXXEventCube())
        test_support._pickeling(self, gammalib.GXXXEventList())
        test_support._pickeling(self, gammalib.GXXXInstDir())
        test_support._pickeling(self, gammalib.GXXXObservation())
        test_support._pickeling(self, gammalib.GXXXResponse())
        test_support._pickeling(self, gammalib.GXXXRoi())

        # Setup test (TODO: to be filled with meaningful values)
        atom = gammalib.GXXXEventAtom()
        bin  = gammalib.GXXXEventBin()
        cube = gammalib.GXXXEventCube()
        list = gammalib.GXXXEventList()
        dir  = gammalib.GXXXInstDir()
        obs  = gammalib.GXXXObservation()
        rsp  = gammalib.GXXXResponse()
        roi  = gammalib.GXXXRoi()

        # Perform pickeling tests of filled classes
        test_support._pickeling(self, gammalib.GXXXEventAtom(atom))
        test_support._pickeling(self, gammalib.GXXXEventBin(bin))
        test_support._pickeling(self, gammalib.GXXXEventCube(cube))
        test_support._pickeling(self, gammalib.GXXXEventList(list))
        test_support._pickeling(self, gammalib.GXXXInstDir(dir))
        test_support._pickeling(self, gammalib.GXXXObservation(obs))
        test_support._pickeling(self, gammalib.GXXXResponse(rsp))
        test_support._pickeling(self, gammalib.GXXXRoi(roi))

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('[INSTRUMENT]')

        # Append tests
        self.append(self._dummy_test, '[INSTRUMENT] dummy test')
        self.append(self._test_pickeling, 'Test [INSTRUMENT] class pickeling')

        # Return
        return
