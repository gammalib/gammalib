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
import os
import sys
import gammalib


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
        Constructor
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Examples directory
        self._cpp_dir    = os.environ['TEST_BUILDDIR'] + '/examples/cpp/'
        self._python_dir = os.environ['TEST_SRCDIR']   + '/examples/python/'

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
        cmd = self._cpp_dir + name

        # Execute binary, make sure we catch any exception
        try:
            rc = os.system(cmd+' > example_'+name+'.log 2>&1')
        except:
            pass

        # Check if execution was successful
        self.test_assert(rc == 0, 'Check "'+name+'" execution on command line')

        # Return
        return

    # Execute Python script
    def _execute_python(self, name):
        """
        Execute Python script
        
        Parameters
        ----------
        name : str
            Python script name with .py extension.
        """
        # Setup command
        cmd = self._python_dir + name + '.py'

        # Execute Python script, make sure we catch any exception
        try:
            rc = os.system(cmd+' > example_'+name+'.log 2>&1')
        except:
            pass

        # Check if execution was successful
        self.test_assert(rc == 0, 'Check "'+name+'" script from command line')

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
        self.append(self.test_matrix_howto, "Test matrix howto")
        self.append(self.test_models_howto, "Test models howto")
        self.append(self.test_xml_howto, "Test XML howto")
        self.append(self.test_xml_html_create, "Test HTML creation")

        # Return
        return

    # Test createfits
    def test_createfits(self):
        """
        Test createfits
        """
        # Execute binary
        self._execute_binary('createfits')

        # Return
        return

    # Test createskymap
    def test_createskymap(self):
        """
        Test createskymap
        """
        # Execute binary
        self._execute_binary('createskymap')

        # Return
        return

    # Test createxml
    def test_createxml(self):
        """
        Test createxml
        """
        # Execute binary
        self._execute_binary('createxml')

        # Return
        return

    # Test interpolate
    def test_interpolate(self):
        """
        Test interpolate
        """
        # Execute binary
        self._execute_binary('interpolate')

        # Return
        return

    # Test numerics
    def test_numerics(self):
        """
        Test numerics
        """
        # Execute binary
        self._execute_binary('numerics')

        # Return
        return

    # Test optimize
    def test_optimize(self):
        """
        Test optimize
        """
        # Execute binary
        self._execute_binary('optimize')

        # Return
        return

    # Test readmodel
    def test_readmodel(self):
        """
        Test readmodel
        """
        # Execute binary
        self._execute_binary('readmodel')

        # Return
        return

    # Test matrix howto
    def test_matrix_howto(self):
        """
        Test matrix howto
        """
        # Execute python script
        self._execute_python('matrix_howto')

        # Return
        return

    # Test models howto
    def test_models_howto(self):
        """
        Test models howto
        """
        # Execute python script
        self._execute_python('models_howto')

        # Return
        return

    # Test XML howto
    def test_xml_howto(self):
        """
        Test XML howto
        """
        # Execute python script
        self._execute_python('xml_howto')

        # Return
        return

    # Test HTML creation
    def test_xml_html_create(self):
        """
        Test HTML creation
        """
        # Execute python script
        self._execute_python('xml_html_create')

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
