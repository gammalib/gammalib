# ==========================================================================
# This module performs unit tests for the GammaLib xspec module.
#
# Copyright (C) 2013 Juergen Knoedlseder
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


# ==================================== #
# Test class for GammaLib xspec module #
# ==================================== #
class Test(GPythonTestSuite):
    """
    Test class for GammaLib xspec module.
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
        self.name("xspec")

        # Append tests
        self.append(self.test_pha, "Test GPha")
        self.append(self.test_arf, "Test GArf")
        self.append(self.test_rmf, "Test GRmf")

        # Return
        return

    # Test GPha class
    def test_pha(self):
        """
        Test GPha class.
        """
        # Allocate GPha object and check its size
        nmeasured = 10
        emeasured = GEbounds(nmeasured, GEnergy(1.0, "TeV"), GEnergy(10.0, "TeV"))
        pha       = GPha(emeasured)
        self.test_value(pha.size(), nmeasured)

        # Set and retrieve a value
        value  = 3.7
        pha[3] = value
        self.test_value(pha[3], value)

        # Return
        return

    # Test GArf class
    def test_arf(self):
        """
        Test GArf class.
        """
        # Allocate GArf object and check its size
        nmeasured = 10
        emeasured = GEbounds(nmeasured, GEnergy(1.0, "TeV"), GEnergy(10.0, "TeV"))
        arf       = GArf(emeasured)
        self.test_value(arf.size(), nmeasured)

        # Set and retrieve a value
        value  = 3.7
        arf[3] = value
        self.test_value(arf[3], value)

        # Return
        return

    # Test GRmf class
    def test_rmf(self):
        """
        Test GRmf class.
        """
        # Allocate GRmf object and check its size
        ntrue     = 20
        nmeasured = 10
        etrue     = GEbounds(ntrue, GEnergy(1.0, "TeV"), GEnergy(10.0, "TeV"))
        emeasured = GEbounds(nmeasured, GEnergy(1.0, "TeV"), GEnergy(10.0, "TeV"))
        rmf       = GRmf(etrue, emeasured)
        self.test_value(rmf.size(), ntrue*nmeasured)
        self.test_value(rmf.ntrue(), ntrue)
        self.test_value(rmf.nmeasured(), nmeasured)

        # Set and retrieve a value
        value    = 3.7
        rmf[3,5] = value
        self.test_value(rmf[3,5], value)
        
        # Return
        return
