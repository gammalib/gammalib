#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the example executables
#
# Copyright (C) 2016 Juergen Knoedlseder
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
import sys
import os


# =============================================== #
# Test class for example executables unit testing #
# =============================================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for GammaLib example executables
    """
    # Constructor
    def __init__(self):
        """
        Constructor.
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Examples directory
        self._dir = '../examples/cpp/'

        # Return
        return

    # Execute binary
    def _execute_binary(self, name):
        """
        Execute binary
        
        Parameters
        ----------
        name : str
            Binary executable name.
        """
        # Setup command
        cmd = self._dir + name

        # Execute binary, make sure we catch any exception
        try:
            rc = os.system(cmd+' >/dev/null 2>&1')
        except:
            pass

        # Check if execution was successful
        self.test_assert(rc == 0, 'Check "'+name+'" execution on command line')

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name("examples")

        # Append tests
        self.append(self.test_createfits, "Test createfits")
        self.append(self.test_createskymap, "Test createskymap")
        self.append(self.test_createxml, "Test createxml")
        self.append(self.test_interpolate, "Test interpolate")
        self.append(self.test_numerics, "Test numerics")
        self.append(self.test_optimize, "Test optimize")
        self.append(self.test_readmodel, "Test readmodel")

        # Return
        return

    # Test createfits
    def test_createfits(self):
        """
        Test createfits
        """
        # Execute binary
        self._execute_binary('createfits')

        #TODO: Do any testing

        # Return
        return

    # Test createskymap
    def test_createskymap(self):
        """
        Test createskymap
        """
        # Execute binary
        self._execute_binary('createskymap')

        #TODO: Do any testing

        # Return
        return

    # Test createxml
    def test_createxml(self):
        """
        Test createxml
        """
        # Execute binary
        self._execute_binary('createxml')

        #TODO: Do any testing

        # Return
        return

    # Test interpolate
    def test_interpolate(self):
        """
        Test interpolate
        """
        # Execute binary
        self._execute_binary('interpolate')

        #TODO: Do any testing

        # Return
        return

    # Test numerics
    def test_numerics(self):
        """
        Test numerics
        """
        # Execute binary
        self._execute_binary('numerics')

        #TODO: Do any testing

        # Return
        return

    # Test optimize
    def test_optimize(self):
        """
        Test optimize
        """
        # Execute binary
        self._execute_binary('optimize')

        #TODO: Do any testing

        # Return
        return

    # Test readmodel
    def test_readmodel(self):
        """
        Test readmodel
        """
        # Execute binary
        self._execute_binary('readmodel')

        #TODO: Do any testing

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Allocate test suites
    suites = gammalib.GTestSuites('Examples testing')

    # Allocate test suite, setup tests and append them to the container
    suite = Test()
    suite.set()
    suites.append(suite)

    # Run test suite
    success = suites.run()

    # Save test results
    suites.save('reports/examples.xml')

    # Set return code
    if success:
        rc = 0
    else:
        rc = 1

    # Exit with return code
    sys.exit(rc)
