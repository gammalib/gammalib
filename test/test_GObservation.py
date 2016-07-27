# ==========================================================================
# This module performs unit tests for the GammaLib observation module.
#
# Copyright (C) 2012-2016 Juergen Knoedlseder
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


# ========================================== #
# Test class for GammaLib observation module #
# ========================================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for GammaLib observation module
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

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('obs')

        # Append tests
        self.append(self.test_energy, 'Test GEnergy class')
        self.append(self.test_energies, 'Test GEnergies class')
        self.append(self.test_time, 'Test GTime class')

        # Return
        return

    # Test GEnergy class
    def test_energy(self):
        """
        Test GEnergy class
        """
        # Setup energies
        energy   = gammalib.GEnergy(3.0, 'TeV')
        energy_b = energy.copy()

        # Unary addition operator
        energy += energy_b
        self.test_value(energy.TeV(), 6.0, 1.0e-7,
             'Check unary addition operator')

        # Unary subtraction operator
        energy -= energy_b
        self.test_value(energy.TeV(), 3.0, 1.0e-7,
             'Check unary subtraction operator')

        # Unary multiplication operator
        energy *= 2.0
        self.test_value(energy.TeV(), 6.0, 1.0e-7,
             'Check unary multiplication operator')

        # Unary division operator
        energy /= 2.0
        self.test_value(energy.TeV(), 3.0, 1.0e-7,
             'Check unary division operator')

        # Binary addition operator
        energy = energy + energy_b
        self.test_value(energy.TeV(), 6.0, 1.0e-7,
             'Check binary addition operator')

        # Binary subtraction operator
        energy = energy - energy_b
        self.test_value(energy.TeV(), 3.0, 1.0e-7,
             'Check binary subtraction operator')

        # Scalar multiplication operator
        energy = energy * 2.0
        self.test_value(energy.TeV(), 6.0, 1.0e-7,
             'Check scalar multiplication operator')

        # Scalar division operator
        energy = energy / 2.0
        self.test_value(energy.TeV(), 3.0, 1.0e-7,
             'Check scalar division operator')

        # Equality operator
        energy1 = gammalib.GEnergy(3.0, 'TeV')
        energy2 = gammalib.GEnergy(5.0, 'TeV')
        self.test_assert(energy1 == energy1, 'Check equality operator')

        # Non-equality operator
        self.test_assert(energy1 != energy2, 'Check non-equality operator')

        # Smaller than operator
        self.test_assert(energy1 < energy2, 'Check smaller than operator')

        # Larger than operator
        self.test_assert(energy2 > energy1, 'Check larger than operator')

        # Smaller or equal than operator
        self.test_assert(energy1 <= energy1, 'Check smaller or equal than operator')
        self.test_assert(energy1 <= energy2, 'Check smaller or equal than operator')

        # Larger or equal than operator
        self.test_assert(energy2 >= energy2, 'Check larger or equal than operator')
        self.test_assert(energy2 >= energy1, 'Check larger or equal than operator')

        # Return
        return

    # Test GEnergies class
    def test_energies(self):
        """
        Test GEnergies class
        """
        # Setup 3 logarithmic energies
        energies = gammalib.GEnergies(3, gammalib.GEnergy(1.0,   'MeV'),
                                         gammalib.GEnergy(100.0, 'MeV'))
        self.test_value(energies[0].MeV(), 1.0, 1.0e-7,
             'Check energy value access')
        self.test_value(energies[1].MeV(), 10.0, 1.0e-7,
             'Check energy value access')
        self.test_value(energies[2].MeV(), 100.0, 1.0e-7,
             'Check energy value access')

        # Return
        return

    # Test GTime class
    def test_time(self):
        """
        Test GTime class
        """
        # Setup times
        time_a = gammalib.GTime(9.0)
        time_b = gammalib.GTime(3.0)
        secs_b = 2.0

        # Unary addition operator
        time  = time_a.copy()
        time += secs_b
        self.test_value(time.secs(), 11.0, 1.0e-7,
             'Check unary addition operator')

        # Unary subtraction operator
        time  = time_a.copy()
        time -= secs_b
        self.test_value(time.secs(), 7.0, 1.0e-7,
             'Check unary subtraction operator')

        # Binary left addition operator
        time = time_a + secs_b
        self.test_value(time.secs(), 11.0, 1.0e-7,
             'Check left addition operator')

        # Binary right addition operator
        time = secs_b + time_a
        self.test_value(time.secs(), 11.0, 1.0e-7,
             'Check right addition operator')

        # Binary subtraction operator
        time = time_a - secs_b
        self.test_value(time.secs(), 7.0, 1.0e-7,
             'Check binary subtraction operator')

        # Binary time subtraction operator
        time = time_a - time_b
        self.test_value(time, 6.0, 1.0e-7,
             'Check binary time subtraction operator')

        # Equality operator
        self.test_assert(time_a == time_a, 'Check equality operator')

        # Non-equality operator
        self.test_assert(time_a != time_b, 'Check non-equality operator')

        # Smaller than operator
        self.test_assert(time_b < time_a, 'Check smaller than operator')

        # Larger than operator
        self.test_assert(time_a > time_b, 'Check larger than operator')

        # Smaller or equal than operator
        self.test_assert(time_b <= time_b,
             'Check smaller or equal than operator')
        self.test_assert(time_b <= time_a,
             'Check smaller or equal than operator')

        # Larger or equal than operator
        self.test_assert(time_a >= time_a,
             'Check larger or equal than operator')
        self.test_assert(time_a >= time_b,
             'Check larger or equal than operator')

        # Return
        return
