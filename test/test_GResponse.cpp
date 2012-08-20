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
void test_lat_response(void)
{
    // Remove FITS file
    system("rm -rf test_rsp.fits");

    std::cout << "Test GLATResponse: Open PSF FITS file: ";
    try {
        // Get HANDOFF Response
        GLATResponse rsp;
        rsp.set_caldb("irf/lat");
        rsp.load("Pass5_v0", "front");

        // Save response
        rsp.save("test_rsp.fits");
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to open FITS file." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ". ok." << std::endl;
}


/***************************************************************************
 *                            Main test function                           *
 ***************************************************************************/
int main(void)
{
    // Dump header
    std::cout << std::endl;
    std::cout << "***************************" << std::endl;
    std::cout << "* GResponse class testing *" << std::endl;
    std::cout << "***************************" << std::endl;

    // Execute the tests
    test_lat_response();

    // Return
    return 0;
}
