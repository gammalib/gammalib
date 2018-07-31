# ==========================================================================
# This module performs unit tests for the GammaLib LAT module.
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
# Test class for GammaLib LAT module #
# ================================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for GammaLib LAT module
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

    # Setup GLATEventList container
    def _setup_eventlist(self):
        """
        Setup GLATEventList container

        Returns
        -------
        list : `~gammalib.GLATEventList`
            GLATEventList container
        """
        # Setup event list container
        list = gammalib.GLATEventList()
        atom = gammalib.GLATEventAtom()
        for i in range(10):
            atom.energy().MeV(float(i))
            list.append(atom)

        # Return event list container
        return list

    # Test GLATEventList class access operators
    def _test_eventlist_access(self):
        """
        Test GLATEventList class observation access
        """
        # Setup event list container
        list = self._setup_eventlist()

        # Perform event list access tests
        test_support.energy_container_access_index(self, list)

        # Return
        return

    # Test GLATEventList class slicing
    def _test_eventlist_slicing(self):
        """
        Test GLATEventList class slicing
        """
        # Setup event list container
        list = self._setup_eventlist()

        # Perform slicing tests
        test_support.energy_container_slicing(self, list)

        # Return
        return

    # Test class pickeling
    def _test_pickeling(self):
        """
        Test class pickeling
        """
        # Set CALDB
        os.environ['CALDB'] = os.environ['TEST_LAT_DATA']+'/../caldb'

        # Perform pickeling tests of empty classes
        test_support.pickeling(self, gammalib.GLATAeff())
        test_support.pickeling(self, gammalib.GLATEdisp())
        #test_support.pickeling(self, gammalib.GLATEventAtom())
        #test_support.pickeling(self, gammalib.GLATEventBin())
        test_support.pickeling(self, gammalib.GLATEventCube())
        test_support.pickeling(self, gammalib.GLATEventList())
        test_support.pickeling(self, gammalib.GLATInstDir())
        test_support.pickeling(self, gammalib.GLATLtCube())
        test_support.pickeling(self, gammalib.GLATObservation())
        test_support.pickeling(self, gammalib.GLATPsf())
        #test_support.pickeling(self, gammalib.GLATResponse())
        test_support.pickeling(self, gammalib.GLATRoi())

        # Setup test
        datadir   = os.environ['TEST_LAT_DATA'] + '/p8v2'
        ltcube    = datadir + '/ltcube.fits'
        eventlist = datadir + '/ft1.fits'
        eventcube = datadir + '/cntmap.fits'
        unbinned  = gammalib.GObservations(datadir + '/obs_unbinned.xml')
        binned    = gammalib.GObservations(datadir + '/obs_binned.xml')
        dir       = gammalib.GSkyDir()
        instdir   = gammalib.GLATInstDir(dir)
        aeff      = os.environ['CALDB']+'/data/glast/lat/bcf/ea/aeff_P8R2_SOURCE_V6_FB.fits'
        edisp     = os.environ['CALDB']+'/data/glast/lat/bcf/edisp/edisp_P8R2_SOURCE_V6_FB.fits'
        psf       = os.environ['CALDB']+'/data/glast/lat/bcf/psf/psf_P8R2_SOURCE_V6_FB.fits'

        # Perform pickeling tests of filled classes
        test_support.pickeling(self, gammalib.GLATAeff(aeff,'FRONT'))
        test_support.pickeling(self, gammalib.GLATEdisp(edisp,'FRONT'))
        #test_support.pickeling(self, gammalib.GLATEventAtom())
        #test_support.pickeling(self, gammalib.GLATEventBin())
        test_support.pickeling(self, gammalib.GLATEventCube(eventcube))
        test_support.pickeling(self, gammalib.GLATEventList(eventlist))
        test_support.pickeling(self, gammalib.GLATInstDir(dir))
        test_support.pickeling(self, gammalib.GLATLtCube(ltcube))
        test_support.pickeling(self, gammalib.GLATObservation(unbinned[0]))
        test_support.pickeling(self, gammalib.GLATObservation(binned[0]))
        test_support.pickeling(self, gammalib.GLATPsf(psf,'FRONT'))
        #test_support.pickeling(self, gammalib.GLATResponse())
        test_support.pickeling(self, gammalib.GLATRoi(instdir, 2.0))

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('LAT')

        # Append tests
        self.append(self._test_eventlist_access, 'Test GLATEventList event access')
        self.append(self._test_eventlist_slicing, 'Test GLATEventList slicing')
        self.append(self._test_pickeling, 'Test LAT class pickeling')

        # Return
        return
