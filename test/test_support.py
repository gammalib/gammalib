# ==========================================================================
# This module provides support functions for the Python unit tests
#
# Copyright (C) 2017-2018 Juergen Knoedlseder
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
import sys
import pickle


# ============================== #
# Test container access by index #
# ============================== #
def _container_access_index(testsuite, container):
    """
    Test container access by index

    Parameters
    ----------
    testsuite : `~gammalib.GPythonTestSuite`
        Test suite class
    container : `~gammalib.GContainer`
        Container class
    """
    # Loop over all elements using the container iterator and count the number
    # of iterations
    niter = 0
    for element in container:
        testsuite.test_value(element.name(), '%s' % niter)
        niter += 1

    # Check that looping was successful
    testsuite.test_value(container.size(), 10, 'Check container size')
    testsuite.test_value(niter, 10, 'Check container iterator')

    # Test access from start
    testsuite.test_value(container[3].name(), '3')

    # Check access from end
    testsuite.test_value(container[-2].name(), '8')

    # Clone first element
    #element = container[0].copy()

    # Check setting by index from start
    #element.name('98')
    #container[3] = element
    #self.test_value(container[3].name(),  '98')

    # Check setting by index from end
    #element.name('99')
    #container[-2] = element
    #self.test_value(container[-2].name(), '99')

    # Return
    return


# ===================================== #
# Test energy container access by index #
# ===================================== #
def _energy_container_access_index(testsuite, container):
    """
    Test energy container access by index

    Parameters
    ----------
    testsuite : `~gammalib.GPythonTestSuite`
        Test suite class
    container : `~gammalib.GContainer`
        Container class
    """
    # Loop over all elements using the container iterator and count the number
    # of iterations
    niter     = 0
    reference = 0.0
    for element in container:
        testsuite.test_value(element.energy().MeV(), reference)
        niter     += 1
        reference += 1.0

    # Check that looping was successful
    testsuite.test_value(container.size(), 10, 'Check container size')
    testsuite.test_value(niter, 10, 'Check container iterator')

    # Test access from start
    testsuite.test_value(container[3].energy().MeV(), 3.0)

    # Check access from end
    testsuite.test_value(container[-2].energy().MeV(), 8.0)

    # Clone first element
    element = container[0].copy()

    # Check setting by index from start
    element.energy().MeV(98.0)
    container[3] = element
    testsuite.test_value(container[3].energy().MeV(),  98.0)

    # Check setting setting by index from end
    element.energy().MeV(99.0)
    container[-2] = element
    testsuite.test_value(container[-2].energy().MeV(),  99.0)

    # Return
    return


# ====================== #
# Test container slicing #
# ====================== #
def _container_slicing(testsuite, container):
    """
    Test container slicing

    Parameters
    ----------
    testsuite : `~gammalib.GPythonTestSuite`
        Test suite class
    container : `~gammalib.GContainer`
        Container class
    """
    # Test container[start:end]
    testsuite.test_value(len(container[3:5]), 2)
    testsuite.test_value(container[3:5][0].name(), '3')
    testsuite.test_value(container[3:5][1].name(), '4')

    # Test container[start:]
    testsuite.test_value(len(container[7:]), 3)
    testsuite.test_value(container[7:][0].name(), '7')
    testsuite.test_value(container[7:][1].name(), '8')
    testsuite.test_value(container[7:][2].name(), '9')

    # Test container[:end]
    testsuite.test_value(len(container[:2]), 2)
    testsuite.test_value(container[:2][0].name(), '0')
    testsuite.test_value(container[:2][1].name(), '1')

    # Test container[:]
    testsuite.test_value(len(container[:]), 10)
    for i in range(10):
        testsuite.test_value(container[:][i].name(), '%s' % i)

    # Test container[start:end:step]
    testsuite.test_value(len(container[3:7:2]), 2)
    testsuite.test_value(container[3:7:2][0].name(), '3')
    testsuite.test_value(container[3:7:2][1].name(), '5')

    # Test container[start:end:step]
    testsuite.test_value(len(container[6:3:-2]), 2)
    testsuite.test_value(container[6:3:-2][0].name(), '6')
    testsuite.test_value(container[6:3:-2][1].name(), '4')

    # Test container[-start:]
    testsuite.test_value(len(container[-2:]), 2)
    testsuite.test_value(container[-2:][0].name(), '8')
    testsuite.test_value(container[-2:][1].name(), '9')

    # Test container[:-end]
    testsuite.test_value(len(container[:-7]), 3)
    testsuite.test_value(container[:-7][0].name(), '0')
    testsuite.test_value(container[:-7][1].name(), '1')
    testsuite.test_value(container[:-7][2].name(), '2')

    # Return
    return


# ============================= #
# Test energy container slicing #
# ============================= #
def _energy_container_slicing(testsuite, container):
    """
    Test energy container slicing

    Parameters
    ----------
    testsuite : `~gammalib.GPythonTestSuite`
        Test suite class
    container : `~gammalib.GContainer`
        Container class
    """
    # Test container[start:end]
    testsuite.test_value(len(container[3:5]), 2)
    testsuite.test_value(container[3:5][0].energy().MeV(), 3.0)
    testsuite.test_value(container[3:5][1].energy().MeV(), 4.0)

    # Test container[start:]
    testsuite.test_value(len(container[7:]), 3)
    testsuite.test_value(container[7:][0].energy().MeV(), 7.0)
    testsuite.test_value(container[7:][1].energy().MeV(), 8.0)
    testsuite.test_value(container[7:][2].energy().MeV(), 9.0)

    # Test container[:end]
    testsuite.test_value(len(container[:2]), 2)
    testsuite.test_value(container[:2][0].energy().MeV(), 0.0)
    testsuite.test_value(container[:2][1].energy().MeV(), 1.0)

    # Test container[:]
    testsuite.test_value(len(container[:]), 10)
    for i in range(10):
        testsuite.test_value(container[:][i].energy().MeV(), float(i))

    # Test container[start:end:step]
    testsuite.test_value(len(container[3:7:2]), 2)
    testsuite.test_value(container[3:7:2][0].energy().MeV(), 3.0)
    testsuite.test_value(container[3:7:2][1].energy().MeV(), 5.0)

    # Test container[start:end:step]
    testsuite.test_value(len(container[6:3:-2]), 2)
    testsuite.test_value(container[6:3:-2][0].energy().MeV(), 6.0)
    testsuite.test_value(container[6:3:-2][1].energy().MeV(), 4.0)

    # Test container[-start:]
    testsuite.test_value(len(container[-2:]), 2)
    testsuite.test_value(container[-2:][0].energy().MeV(), 8.0)
    testsuite.test_value(container[-2:][1].energy().MeV(), 9.0)

    # Test container[:-end]
    testsuite.test_value(len(container[:-7]), 3)
    testsuite.test_value(container[:-7][0].energy().MeV(), 0.0)
    testsuite.test_value(container[:-7][1].energy().MeV(), 1.0)
    testsuite.test_value(container[:-7][2].energy().MeV(), 2.0)

    # Return
    return


# ============== #
# Test pickeling #
# ============== #
def _pickeling(testsuite, object):
    """
    Test class pickeling

    Parameters
    ----------
    testsuite : `~gammalib.GPythonTestSuite`
        Test suite class
    container : `~gammalib.GBase`
        GammaLib class
    """
    # Get class name
    name = object.classname()
    
    # Start exception catching
    testsuite.test_try('Test pickeling of "%s" class.' % name)

    # Test if pickling results in an exception
    try:
        dump       = pickle.dumps(object)
        obj        = pickle.loads(dump)
        dump_again = pickle.dumps(obj)
        testsuite.test_assert(dump == dump_again,
                              'Check pickled dump of "%s\n%s\n%s"' %
                              (name, str(object), str(obj)))
        testsuite.test_try_success()
    except:
        _, e, _ = sys.exc_info()
        testsuite.test_try_failure('Error in pickeling "%s" (%s).' %
                                   (name, str(e)))

    # Return
    return
