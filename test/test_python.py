#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the GammaLib Python bindings.
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
import os
import sys
import gammalib
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# Import test modules
import test_GApplication
import test_GLinalg
import test_GModel
import test_GNumerics
import test_GObservation
import test_GOptimizer
import test_GFits
import test_GSky
import test_GSupport
import test_GTest
import test_GXml
import test_GXspec
import test_GVO

# Try importing MWL tests
try:
    import test_MWL
    has_mwl = True
except:
    has_mwl = False

# Try importing CTA tests
try:
    import test_CTA
    has_cta = True
except:
    has_cta = False

# Try importing LAT tests
try:
    import test_LAT
    has_lat = True
except:
    has_lat = False

# Try importing COMPTEL tests
try:
    import test_COM
    has_com = True
except:
    has_com = False

# Insert new test here

# ================== #
# Perform unit tests #
# ================== #
def test(installed=False):
    """
    Perform unit testing for Python interface
    """
    # Set environment variables first
    if installed:
        os.environ['TEST_DATA']     = 'data'
        os.environ['TEST_COM_DATA'] = 'com/data'
        os.environ['TEST_CTA_DATA'] = 'cta/data'
        os.environ['TEST_LAT_DATA'] = 'lat/data'
        os.environ['TEST_MWL_DATA'] = 'mwl/data'
        # Insert new environment here

    # Allocate test suites
    suites = gammalib.GTestSuites('Python interface testing')

    # Allocate test suites and append them to the container
    suite1  = test_GApplication.Test()
    suite2  = test_GFits.Test()
    suite3  = test_GLinalg.Test()
    suite4  = test_GModel.Test()
    suite5  = test_GNumerics.Test()
    suite6  = test_GObservation.Test()
    suite7  = test_GOptimizer.Test()
    suite8  = test_GSky.Test()
    suite9  = test_GSupport.Test()
    suite10 = test_GTest.Test()
    suite11 = test_GXml.Test()
    suite12 = test_GXspec.Test()
    suite13 = test_GVO.Test()

    # Setup unit tests
    suite1.set()
    suite2.set()
    suite3.set()
    suite4.set()
    suite5.set()
    suite6.set()
    suite7.set()
    suite8.set()
    suite9.set()
    suite10.set()
    suite11.set()
    suite12.set()
    suite13.set()

    # Append tests to container
    suites.append(suite1)
    suites.append(suite2)
    suites.append(suite3)
    suites.append(suite4)
    suites.append(suite5)
    suites.append(suite6)
    suites.append(suite7)
    suites.append(suite8)
    suites.append(suite9)
    suites.append(suite10)
    suites.append(suite11)
    suites.append(suite12)
    suites.append(suite13)

    # Optionally handle MWL suite
    if has_mwl:
        suite_mwl = test_MWL.Test()
        suite_mwl.set()
        suites.append(suite_mwl)

    # Optionally handle CTA suite
    if has_cta:
        suite_cta = test_CTA.Test()
        suite_cta.set()
        suites.append(suite_cta)

    # Optionally handle LAT suite
    if has_lat:
        suite_lat = test_LAT.Test()
        suite_lat.set()
        suites.append(suite_lat)

    # Optionally handle COMPTEL suite
    if has_com:
        suite_com = test_COM.Test()
        suite_com.set()
        suites.append(suite_com)

    # Insert new suite here

    # If we have an installed version then create a temporary
    # directory and copy over all information that is needed
    if installed:

        # Create temporary working directory
        import tempfile
        path = tempfile.mkdtemp()
        os.chdir(path)

        # Get tests directory
        import inspect
        testdir = inspect.getfile(gammalib.tests)
        head, tail = os.path.split(testdir)

        # Copy over test data
        os.system('cp -r %s %s' % (head+'/data', 'data'))
        os.system('cp -r %s %s' % (head+'/cta',  'cta'))

        # Special post processing for CTA files. This is needed because
        # the XML files contain absolute PATH information. This is a kluge
        # that works for now, but it's not a very maintainable way to
        # handle this
        if has_cta:
            xml = gammalib.GXml('cta/data/irf_unbinned.xml')
            elements = xml.element('observation_list').element('observation')
            for element in elements:
                filename   = element.attribute('file')
                head, tail = os.path.split(filename)
                head, dir  = os.path.split(head)
                if len(dir) > 0:
                    filename   = 'cta/'+dir+'/'+tail
                element.attribute('file', filename)
            xml.save('cta/data/irf_unbinned.xml')
            xml = gammalib.GXml('cta/data/irf_1dc.xml')
            elements = xml.element('observation_list').element('observation')
            for element in elements:
                filename   = element.attribute('file')
                head, tail = os.path.split(filename)
                head, dir  = os.path.split(head)
                if dir == 'dc1':
                    dir = 'caldb/dc1'
                if len(dir) > 0:
                    filename = 'cta/'+dir+'/'+tail
                element.attribute('file', filename)
            xml.save('cta/data/irf_1dc.xml')

    # Run test suite
    success = suites.run()

    # If we have a non-installed version then save test results
    if not installed:
        suites.save('reports/python.xml')
    else:
        suites.save('gammalib_reports.xml')

    # Remove temporary direction
    if installed:
        os.system('rm -rf %s' % (path))

    # Raise an exception in case of failure
    if not success:
        raise RuntimeError('At least one error occured during the test.')

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Run tests
    test()
