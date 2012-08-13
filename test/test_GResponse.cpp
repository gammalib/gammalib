/***************************************************************************
 *              test_GResponse.cpp  -  test Response classes               *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#include <stdlib.h>
#include <iostream>
#include "test_GResponse.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */


/***************************************************************************
 *                          Test: LAT Response                             *
 ***************************************************************************/
void test_lat_response(GTestSuite& testsuite)
{
    // Remove FITS file
    system("rm -rf test_rsp.fits");

    testsuite.test_try("Open PSF FITS file");
    try {
        // Get HANDOFF Response
        GLATResponse rsp;
        rsp.set_caldb("irf/lat");
        rsp.load("Pass5_v0", "front");

        // Save response
        rsp.save("test_rsp.fits");

        // All is ok
        testsuit.test_try_success();
    }
    catch (std::exception &e) {
        //There is a problem with the test
        testsuit.test_try_failure(e);
    }
}


/***************************************************************************
 *                            Main test function                           *
 ***************************************************************************/
int main(void)
{
    //Create a test suites container
    GTestSuites testsuites("GResponse class testing");

    // Create a test suite
    GTestSuite testsuite("Test GLATResponse");

    // Append it
    testsuites.append(testsuite);

    // Execute the tests
    test_lat_response(testsuite);

    //End of the tests
    testsuite.test_end();
    // Return
    return 0;
}
