/***************************************************************************
 *                 test_SPI.cpp - Test INTEGRAL/SPI classes                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file test_SPI.cpp
 * @brief INTEGRAL/SPI test class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "test_SPI.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Constants __________________________________________________________ */
const std::string datadir   = std::getenv("TEST_SPI_DATA");
const std::string spi_caldb = datadir + "/../../caldb";
const std::string og_dol    = datadir+"/obs/og_spi.fits";


/***********************************************************************//**
 * @brief Set INTEGRAL/SPI test methods
 ***************************************************************************/
void TestGSPI::set(void)
{
    // Set test name
    name("GSPI");

    // Append tests to test suite
    append(static_cast<pfunction>(&TestGSPI::test_obs),
           "Test GSPIObservation");
    append(static_cast<pfunction>(&TestGSPI::test_response),
           "Test GSPIResponse");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGSPI* TestGSPI::clone(void) const
{
    // Clone test suite
    return new TestGSPI(*this);
}


/***********************************************************************//**
 * @brief Test GSPIObservation class
 ***************************************************************************/
void TestGSPI::test_obs(void)
{
    // Test empty GSPIObservation instance
    GSPIObservation obs1;

    // Construct GSPIObservation instance from Observation group
    GSPIObservation obs2(og_dol);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GSPIResponse class
 ***************************************************************************/
void TestGSPI::test_response(void)
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
    GTestSuites testsuites("INTEGRAL/SPI instrument specific class testing");

    // Set CALDB environment variable
    std::string caldb = "CALDB="+spi_caldb;
    putenv((char*)caldb.c_str());

    // Initially assume that we pass all tests
    bool success = true;

    // Create test suites and append them to the container
    TestGSPI suite;
    testsuites.append(suite);

    // Run the testsuites
    success = testsuites.run();

    // Save test report
    testsuites.save("reports/GSPI.xml");

    // Return success status
    return (success ? 0 : 1);
}
