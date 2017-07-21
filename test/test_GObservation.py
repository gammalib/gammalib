# ==========================================================================
# This module performs unit tests for the GammaLib observation module.
#
# Copyright (C) 2012-2017 Juergen Knoedlseder
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

    # Setup GObservations container
    def _setup_obs(self):
        """
        Setup GObservations container

        Returns
        -------
        obs : `~gammalib.GObservations`
            Observation container
        """
        # Setup observation container
        obs = gammalib.GObservations()
        run = gammalib.GCTAObservation()
        for i in range(10):
            run.name('%s' % i)
            run.id('%s' % i)
            obs.append(run)

        # Setup model container
        models = gammalib.GModels()
        model  = gammalib.GModelSky()
        for i in range(10):
            model.name('%s' % i)
            models.append(model)

        # Set observation's model container
        obs.models(models)

        # Return observations container
        return obs

    # Test GEnergy class
    def _test_energy(self):
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
    def _test_energies(self):
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
    def _test_time(self):
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

    # Test GObservations class access operators
    def _test_observations_access(self):
        """
        Test GObservations class observation access
        """
        # Setup observation container and CTA observation
        obs = self._setup_obs()
        run = gammalib.GCTAObservation()

        # Perform observation access tests
        test_support._container_access_index(self, obs)

        # Check observation setting by index from start
        run.id('98')
        obs[3] = run
        self.test_value(obs[3].id(),  '98')

        # Check observation setting by index from end
        run.id('99')
        obs[-2] = run
        self.test_value(obs[-2].id(), '99')

        # Return
        return

    # Test GObservations class slicing
    def _test_observations_slicing(self):
        """
        Test GObservations class slicing
        """
        # Setup observation container
        obs = self._setup_obs()

        # Perform slicing tests
        test_support._container_slicing(self, obs)

        # Perform additional slicing tests
        self.test_value(len(obs[3:5].models()), 10)
        self.test_value(len(obs[7:].models()), 10)
        self.test_value(len(obs[:2].models()), 10)
        self.test_value(len(obs[:].models()), 10)
        self.test_value(len(obs[3:7:2].models()), 10)
        self.test_value(len(obs[6:3:-2].models()), 10)
        self.test_value(len(obs[-2:].models()), 10)
        self.test_value(len(obs[:-7].models()), 10)

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
        self.append(self._test_energy, 'Test GEnergy class')
        self.append(self._test_energies, 'Test GEnergies class')
        self.append(self._test_time, 'Test GTime class')
        self.append(self._test_observations_access, 'Test GObservations class observation access')
        self.append(self._test_observations_slicing, 'Test GObservations class slicing')

        # Return
        return
