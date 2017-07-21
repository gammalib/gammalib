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

    # Setup GTimes container
    def _setup_times(self):
        """
        Setup GTimes container

        Returns
        -------
        times : `~gammalib.GTimes`
            GTime container
        """
        # Setup time container
        times = gammalib.GTimes()
        time  = gammalib.GTime()
        for i in range(10):
            time.secs(float(i))
            times.append(time)

        # Return times container
        return times

    # Setup GPhotons container
    def _setup_photons(self):
        """
        Setup GPhotons container

        Returns
        -------
        photons : `~gammalib.GPhotons`
            GPhoton container
        """
        # Setup photon container
        photons = gammalib.GPhotons()
        photon  = gammalib.GPhoton()
        for i in range(10):
            photon.energy().MeV(float(i))
            photons.append(photon)

        # Return photons container
        return photons

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
        # Setup energies container
        energies = gammalib.GEnergies(3, gammalib.GEnergy(1.0,   'MeV'),
                                         gammalib.GEnergy(100.0, 'MeV'))

        # Check energies container
        self.test_value(energies[0].MeV(), 1.0, 1.0e-7,
             'Check energy value access')
        self.test_value(energies[1].MeV(), 10.0, 1.0e-7,
             'Check energy value access')
        self.test_value(energies[2].MeV(), 100.0, 1.0e-7,
             'Check energy value access')

        # Return
        return

    # Test GEnergies class access operators
    def _test_energies_access(self):
        """
        Test GEnergies class observation access
        """
        # Setup energies container
        energies = gammalib.GEnergies(10, gammalib.GEnergy(1.0,  'MeV'),
                                          gammalib.GEnergy(10.0, 'MeV'), False)

        # Loop over all elements using the container iterator and count the
        # number of iterations
        nenergies = 0
        reference = 1.0
        for energy in energies:
            self.test_value(energy.MeV(), reference)
            nenergies += 1
            reference += 1.0

        # Check that looping was successful
        self.test_value(energies.size(), 10, 'Check container size')
        self.test_value(nenergies, 10, 'Check container iterator')

        # Test access from start
        self.test_value(energies[3].MeV(), 4.0)

        # Check access from end
        self.test_value(energies[-2].MeV(), 9.0)

        # Return
        return

    # Test GEnergies class slicing
    def _test_energies_slicing(self):
        """
        Test GEnergies class slicing
        """
        # Setup energies container
        energies = gammalib.GEnergies(10, gammalib.GEnergy(1.0,  'MeV'),
                                          gammalib.GEnergy(10.0, 'MeV'), False)

        # Test energies[start:end]
        self.test_value(len(energies[3:5]), 2)
        self.test_value(energies[3:5][0].MeV(), 4.0)
        self.test_value(energies[3:5][1].MeV(), 5.0)

        # Test energies[start:]
        self.test_value(len(energies[7:]), 3)
        self.test_value(energies[7:][0].MeV(), 8.0)
        self.test_value(energies[7:][1].MeV(), 9.0)
        self.test_value(energies[7:][2].MeV(), 10.0)

        # Test energies[:end]
        self.test_value(len(energies[:2]), 2)
        self.test_value(energies[:2][0].MeV(), 1.0)
        self.test_value(energies[:2][1].MeV(), 2.0)

        # Test energies[:]
        self.test_value(len(energies[:]), 10)
        for i in range(10):
            self.test_value(energies[:][i].MeV(), float(i+1))

        # Test energies[start:end:step]
        self.test_value(len(energies[3:7:2]), 2)
        self.test_value(energies[3:7:2][0].MeV(), 4.0)
        self.test_value(energies[3:7:2][1].MeV(), 6.0)

        # Test energies[start:end:step]
        self.test_value(len(energies[6:3:-2]), 2)
        self.test_value(energies[6:3:-2][0].MeV(), 7.0)
        self.test_value(energies[6:3:-2][1].MeV(), 5.0)

        # Test energies[-start:]
        self.test_value(len(energies[-2:]), 2)
        self.test_value(energies[-2:][0].MeV(), 9.0)
        self.test_value(energies[-2:][1].MeV(), 10.0)

        # Test energies[:-end]
        self.test_value(len(energies[:-7]), 3)
        self.test_value(energies[:-7][0].MeV(), 1.0)
        self.test_value(energies[:-7][1].MeV(), 2.0)
        self.test_value(energies[:-7][2].MeV(), 3.0)

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

    # Test GTimes class access operators
    def _test_times_access(self):
        """
        Test GTimes class observation access
        """
        # Setup times container
        times = self._setup_times()

        # Loop over all elements using the container iterator and count the
        # number of iterations
        ntimes    = 0
        reference = 0.0
        for time in times:
            self.test_value(time.secs(), reference)
            ntimes    += 1
            reference += 1.0

        # Check that looping was successful
        self.test_value(times.size(), 10, 'Check container size')
        self.test_value(ntimes, 10, 'Check container iterator')

        # Test access from start
        self.test_value(times[3].secs(), 3.0)

        # Check access from end
        self.test_value(times[-2].secs(), 8.0)

        # Return
        return

    # Test GTimes class slicing
    def _test_times_slicing(self):
        """
        Test GTimes class slicing
        """
        # Setup times container
        times = self._setup_times()

        # Test times[start:end]
        self.test_value(len(times[3:5]), 2)
        self.test_value(times[3:5][0].secs(), 3.0)
        self.test_value(times[3:5][1].secs(), 4.0)

        # Test times[start:]
        self.test_value(len(times[7:]), 3)
        self.test_value(times[7:][0].secs(), 7.0)
        self.test_value(times[7:][1].secs(), 8.0)
        self.test_value(times[7:][2].secs(), 9.0)

        # Test times[:end]
        self.test_value(len(times[:2]), 2)
        self.test_value(times[:2][0].secs(), 0.0)
        self.test_value(times[:2][1].secs(), 1.0)

        # Test times[:]
        self.test_value(len(times[:]), 10)
        for i in range(10):
            self.test_value(times[:][i].secs(), float(i))

        # Test times[start:end:step]
        self.test_value(len(times[3:7:2]), 2)
        self.test_value(times[3:7:2][0].secs(), 3.0)
        self.test_value(times[3:7:2][1].secs(), 5.0)

        # Test times[start:end:step]
        self.test_value(len(times[6:3:-2]), 2)
        self.test_value(times[6:3:-2][0].secs(), 6.0)
        self.test_value(times[6:3:-2][1].secs(), 4.0)

        # Test times[-start:]
        self.test_value(len(times[-2:]), 2)
        self.test_value(times[-2:][0].secs(), 8.0)
        self.test_value(times[-2:][1].secs(), 9.0)

        # Test times[:-end]
        self.test_value(len(times[:-7]), 3)
        self.test_value(times[:-7][0].secs(), 0.0)
        self.test_value(times[:-7][1].secs(), 1.0)
        self.test_value(times[:-7][2].secs(), 2.0)

        # Return
        return

    # Test GPhotons class access operators
    def _test_photons_access(self):
        """
        Test GPhotons class observation access
        """
        # Setup photons container
        photons = self._setup_photons()

        # Loop over all elements using the container iterator and count the
        # number of iterations
        nphotons  = 0
        reference = 0.0
        for photon in photons:
            self.test_value(photon.energy().MeV(), reference)
            nphotons  += 1
            reference += 1.0

        # Check that looping was successful
        self.test_value(photons.size(), 10, 'Check container size')
        self.test_value(nphotons, 10, 'Check container iterator')

        # Test access from start
        self.test_value(photons[3].energy().MeV(), 3.0)

        # Check access from end
        self.test_value(photons[-2].energy().MeV(), 8.0)

        # Return
        return

    # Test GPhotons class slicing
    def _test_photons_slicing(self):
        """
        Test GPhotons class slicing
        """
        # Setup photons container
        photons = self._setup_photons()

        # Test photons[start:end]
        self.test_value(len(photons[3:5]), 2)
        self.test_value(photons[3:5][0].energy().MeV(), 3.0)
        self.test_value(photons[3:5][1].energy().MeV(), 4.0)

        # Test photons[start:]
        self.test_value(len(photons[7:]), 3)
        self.test_value(photons[7:][0].energy().MeV(), 7.0)
        self.test_value(photons[7:][1].energy().MeV(), 8.0)
        self.test_value(photons[7:][2].energy().MeV(), 9.0)

        # Test photons[:end]
        self.test_value(len(photons[:2]), 2)
        self.test_value(photons[:2][0].energy().MeV(), 0.0)
        self.test_value(photons[:2][1].energy().MeV(), 1.0)

        # Test photons[:]
        self.test_value(len(photons[:]), 10)
        for i in range(10):
            self.test_value(photons[:][i].energy().MeV(), float(i))

        # Test photons[start:end:step]
        self.test_value(len(photons[3:7:2]), 2)
        self.test_value(photons[3:7:2][0].energy().MeV(), 3.0)
        self.test_value(photons[3:7:2][1].energy().MeV(), 5.0)

        # Test photons[start:end:step]
        self.test_value(len(photons[6:3:-2]), 2)
        self.test_value(photons[6:3:-2][0].energy().MeV(), 6.0)
        self.test_value(photons[6:3:-2][1].energy().MeV(), 4.0)

        # Test photons[-start:]
        self.test_value(len(photons[-2:]), 2)
        self.test_value(photons[-2:][0].energy().MeV(), 8.0)
        self.test_value(photons[-2:][1].energy().MeV(), 9.0)

        # Test photons[:-end]
        self.test_value(len(photons[:-7]), 3)
        self.test_value(photons[:-7][0].energy().MeV(), 0.0)
        self.test_value(photons[:-7][1].energy().MeV(), 1.0)
        self.test_value(photons[:-7][2].energy().MeV(), 2.0)

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
        self.append(self._test_energy, 'Test GEnergy')
        self.append(self._test_energies, 'Test GEnergies')
        self.append(self._test_energies_access, 'Test GEnergies energy access')
        self.append(self._test_energies_slicing, 'Test GEnergies slicing')
        self.append(self._test_time, 'Test GTime')
        self.append(self._test_times_access, 'Test GTimes time access')
        self.append(self._test_times_slicing, 'Test GTimes slicing')
        self.append(self._test_photons_access, 'Test GPhotons photon access')
        self.append(self._test_photons_slicing, 'Test GPhotons slicing')
        self.append(self._test_observations_access, 'Test GObservations observation access')
        self.append(self._test_observations_slicing, 'Test GObservations slicing')

        # Return
        return
