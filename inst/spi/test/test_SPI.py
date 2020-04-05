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
import os
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

        # Set test directories
        self._data  = os.environ['TEST_SPI_DATA']
        self._caldb = self._data + '/../caldb'

        # Return
        return

    # Test GSPIEventCube class
    def _event_cube_test(self):
        """
        Test GSPIEventCube class
        """
        # Set Observation Group file name
        og_dol = self._data + '/obs/og_spi.fits'

        # Load event cube
        cube = gammalib.GSPIEventCube(og_dol)
        self.test_value(cube.size(), 16720, 'Check cube size')
        self.test_value(cube.dim(), 3, 'Check cube dimension')
        self.test_value(cube.naxis(0), 88, 'Check number of pointings')
        self.test_value(cube.naxis(1), 19, 'Check number of detectors')
        self.test_value(cube.naxis(2), 10, 'Check number of energy bins')
        self.test_value(cube.number(), 101269457, 'Check number of events')

        # Determine total number of counts and model events
        counts = 0.0
        model  = 0.0
        for event in cube:
            counts += event.counts()
            model  += event.model(0)
        self.test_value(counts, 101269457.0,      1.0e-6, 'Check number of counts')
        self.test_value(model,  101269456.964783, 1.0e-6, 'Check model counts')

        # Return
        return

    # Test GSPIObservation class
    def _observation_test(self):
        """
        Test GSPIObservation class
        """
        # Set Observation Group file name
        og_dol = self._data + '/obs/og_spi.fits'

        # Load event cube
        obs = gammalib.GSPIObservation(og_dol)
        self.test_value(obs.ontime(), 193966.8178673, 1.0e-6, 'Check ontime')
        self.test_value(obs.livetime(), 170657.5371606, 1.0e-6, 'Check livetime')
        self.test_value(obs.deadc(), 0.8798285, 1.0e-6, 'Check deadtime correction')

        # Set XML file name
        xml = self._data + '/obs.xml'

        # Load observation from XML file
        obs = gammalib.GObservations(xml)
        self.test_value(obs[0].name(), 'Crab', 'Check observation name')
        self.test_value(obs[0].id(), '0044', 'Check observation identifier')
        self.test_value(obs[0].instrument(), 'SPI', 'Check instrument name')
        self.test_value(obs[0].ontime(), 193966.8178673, 1.0e-6, 'Check ontime')
        self.test_value(obs[0].livetime(), 170657.5371606, 1.0e-6, 'Check livetime')
        self.test_value(obs[0].deadc(), 0.8798285, 1.0e-6, 'Check deadtime correction')

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
        self.append(self._event_cube_test, 'Test GSPIEventCube class')
        self.append(self._observation_test, 'Test GSPIObservation class')
        self.append(self._test_pickeling, 'Test INTEGRAL/SPI class pickeling')

        # Return
        return
