/***************************************************************************
 *                   test_GVO.hpp - Test VO module                         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
 * @file test_GVO.hpp
 * @brief Implementation of unit tests for VO module
 * @author Juergen Knoedlseder 
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "test_GVO.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @brief Set parameters and tests
 **************************************************************************/
void TestGVO::set(void)
{
    // Test name
    name("GVO");

    // Append tests
    append(static_cast<pfunction>(&TestGVO::test_GVOClient), "Test GVOClient class");

    // Return
    return; 
}


/***********************************************************************//**
 * @brief Test VO client class
 **************************************************************************/
void TestGVO::test_GVOClient(void)
{
    // Test constructor
    test_try("GVOClient empty constructor");
    try {
        GVOClient client;
        std::cout << client << std::endl;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test Hub connection
    test_try("GVOClient hub connection and disconnection");
    try {
        GVOClient client;
        client.connect();
        std::cout << client << std::endl;
        client.disconnect();
        std::cout << client << std::endl;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Return
    return;
}


/***************************************************************************
 * @brief Main entry point for test executable
 ***************************************************************************/
int main(void)
{
    // Allocate test suit container
    GTestSuites testsuites("VO module");

    // Initially assume that we pass all tests
    bool success = true;

    // Create a test suite
    TestGVO test;

    // Append test suite to the container
    testsuites.append(test);

    // Run the testsuites
    success = testsuites.run();

    // Save test report
    testsuites.save("reports/GVO.xml");

    // Return success status
    return (success ? 0 : 1);
}
