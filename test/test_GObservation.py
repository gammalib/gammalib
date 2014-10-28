# ==========================================================================
# This module performs unit tests for the GammaLib observation module.
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


# ========================================== #
# Test class for GammaLib observation module #
# ========================================== #
class Test(GPythonTestSuite):
    """
    Test class for GammaLib observation module.
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
        self.name("obs")

        # Append tests
        self.append(self.test_energy, "Test GEnergy class")
        self.append(self.test_time, "Test GTime class")

        # Return
        return

    # Test GEnergy class
    def test_energy(self):
        """
        Test GEnergy class.
        """
        # Setup energies
        energy    = GEnergy(3.0, "TeV")
        energy_b  = energy.copy()

        # Unary addition operator
        energy += energy_b
        self.test_value(energy.TeV(), 6.0);

        # Unary subtraction operator
        energy -= energy_b
        self.test_value(energy.TeV(), 3.0);

        # Binary addition operator
        energy = energy + energy_b
        self.test_value(energy.TeV(), 6.0);

        # Binary subtraction operator
        energy = energy - energy_b
        self.test_value(energy.TeV(), 3.0);

        # Scalar multiplication operator
        energy = energy * 2.0
        self.test_value(energy.TeV(), 6.0);

        # Scalar division operator
        energy = energy / 2.0
        self.test_value(energy.TeV(), 3.0);

        # Return
        return


    # Test GTime class
    def test_time(self):
        """
        Test GTime class.
        """
        # Setup times
        time   = GTime(3.0)
        time_b = time.copy()

        # Binary addition operator
        time = time + time_b
        self.test_value(time.secs(), 6.0);

        # Binary subtraction operator
        time = time - time_b
        self.test_value(time.secs(), 3.0);

        # Scalar multiplication operator
        time = time * 2.0
        self.test_value(time.secs(), 6.0);

        # Scalar division operator
        time = time / 2.0
        self.test_value(time.secs(), 3.0);

        # Return
        return
