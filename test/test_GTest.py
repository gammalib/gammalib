# ==========================================================================
# This module performs unit tests for the GammaLib test module.
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
import gammalib
import test_support


# =================================== #
# Test class for GammaLib test module #
# =================================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for GammaLib test module
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

    # Setup test suite container
    def _setup_suites(self):
        """
        Setup test suite container

        Returns
        -------
        suites : `~gammalib.GTestSuites`
            Test suite container
        """
        # Setup test suite container
        suites = gammalib.GTestSuites()
        suite  = gammalib.GPythonTestSuite()
        for i in range(10):
            suite.name('%s' % i)
            suites.append(suite)

        # Return test suite container
        return suites

    # Test GTestSuites class access operators
    def _test_suites_access(self):
        """
        Test GTestSuites class parameter access
        """
        # Setup test suite container and test suite
        suites = self._setup_suites()
        suite  = gammalib.GPythonTestSuite()

        # Perform GTestSuites access tests
        test_support.container_access_index(self, suites)

        # Check parameter setting by index from start
        suite.name('98')
        suites[3] = suite
        self.test_value(suites[3].name(), '98')

        # Check parameter setting by index from end
        suite.name('99')
        suites[-2] = suite
        self.test_value(suites[-2].name(), '99')

        # Return
        return

    # Test GTestSuites class slicing
    def _test_suites_slicing(self):
        """
        Test GTestSuites class slicing
        """
        # Setup test suite container
        suites = self._setup_suites()

        # Perform slicing tests
        test_support.container_slicing(self, suites)

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('test')

        # Append tests
        self.append(self._test_suites_access, 'Test GTestSuites parameter access')
        self.append(self._test_suites_slicing, 'Test GTestSuites slicing')

        # Return
        return
