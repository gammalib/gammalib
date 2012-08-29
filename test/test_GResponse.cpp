/***************************************************************************
 *              test_GResponse.cpp  -  test Response classes               *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008-2012 by Jurgen Knodlseder              *
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

void TestGResponse::set(void){
    // Test name
    name("GReponse);

    //add tests
    add_test(static_cast<pfunction>(&TestGResponse::test_lat_response),"Test LAT Response");

    return;
}
/***************************************************************************
 *                          Test: LAT Response                             *
 ***************************************************************************/
void TestGResponse::test_lat_response(void)
{
    // Remove FITS file
    system("rm -rf test_rsp.fits");

    // Open PSF FITS file

    // Get HANDOFF Response
    GLATResponse rsp;
    rsp.set_caldb("irf/lat");
    rsp.load("Pass5_v0", "front");

    // Save response
    rsp.save("test_rsp.fits");

    return;
}


/***************************************************************************
 *                            Main test function                           *
 ***************************************************************************/
int main(void)
{
    GTestSuites testsuites("GResponse");

    bool was_successful=true;

    //Create a test suite
    TestGReponse test;

    //Append to the container
    testsuites.append(test);

    //Run
    was_successful=testsuites.run();

    //save xml report
    testsuites.save("reports/GResponse.xml");

    // Return
    return was_successful ? 0:1;
}
