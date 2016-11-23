#! /usr/bin/env python
# ==========================================================================
# Display test reports
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
import glob
import gammalib


# ====================== #
# Report test case error #
# ====================== #
def report_one_error(name, xml):
    """
    Report test case error

    Parameters
    ----------
    name : str
        Test case name
    xml : gammalib.GXmlElement
        Error XML element
    """
    # Extract information
    message = xml.attribute('message')
    type    = xml.attribute('type')

    # Print error information
    print(' Error in %s: %s (%s)' % (name, message, type))

    # Return
    return


# ======================== #
# Report test case failure #
# ======================== #
def report_one_failure(name, xml):
    """
    Report test case failure

    Parameters
    ----------
    name : str
        Test case name
    xml : gammalib.GXmlElement
        Failure XML element
    """
    # Extract information
    message = xml.attribute('message')
    type    = xml.attribute('type')

    # Print failure information
    print(' Failure in %s: %s (%s)' % (name, message, type))

    # Return
    return


# ======================= #
# Report test case result #
# ======================= #
def report_one_test_case(xml):
    """
    Report test case result

    Parameters
    ----------
    xml : gammalib.GXmlElement
        Test case XML element
    """
    # Extract test case information
    name     = xml.attribute('name')
    errors   = int(xml.elements('error'))
    failures = int(xml.elements('failure'))

    # If there are errors then display them
    if errors > 0:
        num = int(xml.elements('error'))
        for i in range(num):
            error = xml.element('error[%d]' % i)
            report_one_error(name, error)

    # If there are failures then display them
    if failures > 0:
        num = int(xml.elements('failure'))
        for i in range(num):
            failure = xml.element('failure[%d]' % i)
            report_one_failure(name, failure)

    # Return
    return


# ========================= #
# Report test suite results #
# ========================= #
def report_one_test_suite(xml):
    """
    Report test suite result

    Parameters
    ----------
    xml : gammalib.GXmlElement
        Test suite XML element
    """
    # Extract test suite information
    name     = xml.attribute('name')
    errors   = int(xml.attribute('errors'))
    failures = int(xml.attribute('failures'))

    # If there are errors or failures then print them
    if errors > 0 or failures > 0:

        # Print header
        if errors == 0:
            print('%d failures occured in module "%s":' %
                  (failures, name))
        elif failures == 0:
            print('%d errors occured in module "%s":' %
                  (errors, name))
        else:
            print('%d errors and %d failures occured in module "%s":' %
                  (errors, failures, name))

        # Determine number of test cases
        num = int(xml.elements('testcase'))

        # Loop over all test cases
        for i in range(num):

            # Get test case
            case = xml.element('testcase[%d]' % i)

            # Report test case result
            report_one_test_case(case)

    # Return
    return


# ====================== #
# Report one test result #
# ====================== #
def report_one_test_result(file):
    """
    Report one test result

    Parameters
    ----------
    file : str
        Test result file
    """
    # Open XML file
    xml = gammalib.GXml(file)

    # Determine number of testsuites
    num = int(xml.element('testsuites').elements('testsuite'))

    # Loop over all test suites
    for i in range(num):

        # Get test suite
        suite = xml.element('testsuites > testsuite[%d]' % i)

        # Report test suite result
        report_one_test_suite(suite)

    # Return
    return


# =================== #
# Report test results #
# =================== #
def report_test_results(path, pattern="*.xml"):
    """
    Report test results

    Parameters
    ----------
    path : str
        Directory containing test results
    pattern : str, optional
        Test report pattern
    """
    # Get list of test reports
    reports = glob.glob(path+'/'+pattern)

    # Loop over all reports
    for report in reports:
        report_one_test_result(report)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Print header
    print('')
    print('GammaLib test reports')
    print('=====================')
    
    # Report test results
    report_test_results('test/reports')
