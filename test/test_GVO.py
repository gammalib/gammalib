# ==========================================================================
# This module performs unit tests for the GammaLib vo module.
#
# Copyright (C) 2014 Juergen Knoedlseder
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


# ================================= #
# Test class for GammaLib vo module #
# ================================= #
class Test(GPythonTestSuite):
    """
    Test class for GammaLib vo module.
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
        self.name("vo")

        # Append tests
        self.append(self.test_voclient, "Test GVOClient")

        # Return
        return

    # Test GVOClient class
    def test_voclient(self):
        """
        Test GVOClient class.
        """
        # Allocate GVOClient object
        client = GVOClient()

        # Return
        return
