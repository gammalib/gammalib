/***************************************************************************
 *           test_GApplication.cpp  -  test GApplication classes           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
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
/**
 * @file test_GApplication.cpp
 * @brief Testing of application classes.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GammaLib.hpp"
#include <iostream>
#include <stdexcept>
#include <stdlib.h>

/* __ Globals ____________________________________________________________ */


/***************************************************************************
 * Test: GLog                                                              *
 ***************************************************************************/
void test_GLog(GTestSuite& testsuite)
{
    // Test constructor
    testsuite.test_try("Test GLog constructor");
    try {
        GLog log1;
        GLog log2 = log1;

        // Test success
        testsuite.test_try_success();
    }
    catch (std::exception &e) {
        testsuite.test_try_failure(e);
    }

    // Test stream logger
    testsuite.test_try("Test stream logger");
    try {
        GLog log;
        log.date(true);
        log.name("Test");
        log.open("test_GApplication.log", true);
        log << "1. This is a C++ string" << std::endl;
        log << "2. This is an integer: " << int(1) << std::endl;
        log << "3. This is a single precision: " << float(1.23456789) << std::endl;
        log << "4. This is a double precision: " << double(1.23456789) << std::endl;
        log << "5. This is a character: " << 'a' << std::endl;
        log << "6. This is a C string: " << "a" << std::endl;
        log << "7. This is a Boolean: " << true << std::endl;

        testsuite.test_try_success();
    }
    catch (std::exception &e) {
        testsuite.test_try_failure(e);
    }

    // Test C logger
    testsuite.test_try("Test C logger");
    try {
        GLog log;
        log.date(true);
        log.name("Test");
        log.open("test_GApplication.log");
        log("%s", "8. This is a C++ string");
        log("%s %d", "9. This is an integer:", int(1));
        log("%s %f", "10. This is a single precision:", float(1.23456789));
        log("%s %f", "11. This is a double precision:", double(1.23456789));

        testsuite.test_try_success();
    }
    catch (std::exception &e) {
        testsuite.test_try_failure(e);
    }

    // Exit test
    return;

}


/***************************************************************************
 *                            Main test function                           *
 ***************************************************************************/
int main(void)
{
    // Create a test suite container
    GTestSuites testsuites("GApplication classes testing");

    // Create a test suite
    GTestSuite testsuite("Test GLog");

    // Append testsuite
    testsuites.append(testsuite);

    // Execute tests
    test_GLog(testsuite);

    // This test suite is finished
    testsuite.end_test();

    //save xml report
    testsuites.save("reports/GApplication.xml");

    // Return
    return 0;
}
