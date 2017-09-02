/***************************************************************************
 *                 test_XXX.cpp - Test [INSTRUMENT] classes                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) [YEAR] by [AUTHOR]                                       *
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
 * @file test_XXX.cpp
 * @brief [INSTRUMENT] test class implementation
 * @author [AUTHOR]
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "test_XXX.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Constants __________________________________________________________ */
const std::string datadir   = std::getenv("TEST_XXX_DATA");
const std::string xxx_caldb = datadir + "/../../caldb";


/***********************************************************************//**
 * @brief Set [INSTRUMENT] test methods
 ***************************************************************************/
void TestGXXX::set(void)
{
    // Set test name
    name("GXXX");

    // Append tests to test suite
    append(static_cast<pfunction>(&TestGXXX::test_obs),
           "Test GXXXObservation");
    append(static_cast<pfunction>(&TestGXXX::test_response),
           "Test GXXXResponse");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGXXX* TestGXXX::clone(void) const
{
    // Clone test suite
    return new TestGXXX(*this);
}


/***********************************************************************//**
 * @brief Test GXXXObservation class
 ***************************************************************************/
void TestGXXX::test_obs(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GXXXResponse class
 ***************************************************************************/
void TestGXXX::test_response(void)
{
    // Return
    return;
}


/***************************************************************************
 * @brief Main entry point for test executable
 ***************************************************************************/
int main(void)
{
    // Allocate test suit container
    GTestSuites testsuites("[INSTRUMENT] instrument specific class testing");

    // Set CALDB environment variable
    std::string caldb = "CALDB="+xxx_caldb;
    putenv((char*)caldb.c_str());

    // Initially assume that we pass all tests
    bool success = true;

    // Create test suites and append them to the container
    TestGXXX suite;
    testsuites.append(suite);

    // Run the testsuites
    success = testsuites.run();

    // Save test report
    testsuites.save("reports/GXXX.xml");

    // Return success status
    return (success ? 0 : 1);
}
