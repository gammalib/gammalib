/***************************************************************************
 *                      test_SPI.cpp  -  test SPI classes                  *
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
 * @file test_SPI.cpp
 * @brief Testing of SPI classes
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <unistd.h>
#include "GTools.hpp"
#include "test_SPI.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Constants __________________________________________________________ */
const std::string datadir       = "../inst/com/test/data";
const std::string com_caldb     = "../inst/com/caldb";
const std::string com_iaq       = "u47569_iaq.fits";          // 1-3 MeV
const std::string com_dre       = datadir+"/m50439_dre.fits"; // 1-3 MeV
const std::string com_drb       = datadir+"/m34997_drg.fits";
const std::string com_drg       = datadir+"/m34997_drg.fits";
const std::string com_drx       = datadir+"/m32171_drx.fits";
const std::string com_obs       = datadir+"/obs.xml";
const std::string com_model     = datadir+"/crab.xml";


/***********************************************************************//**
 * @brief Set SPI response test methods
 ***************************************************************************/
void TestGSPIResponse::set(void)
{
    // Set test name
    name("GSPIResponse");

    // Append tests to test suite
    append(static_cast<pfunction>(&TestGSPIResponse::test_inst_dir), "Test instrument direction");
    append(static_cast<pfunction>(&TestGSPIResponse::test_pointing), "Test pointing");
    append(static_cast<pfunction>(&TestGSPIResponse::test_response), "Test response");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set SPI observation test methods
 ***************************************************************************/
void TestGSPIObservation::set(void)
{
    // Set test name
    name("GSPIObservation");

    // Append tests to test suite
    append(static_cast<pfunction>(&TestGSPIObservation::test_binned_obs), "Test binned observation");
    append(static_cast<pfunction>(&TestGSPIObservation::test_event_bin), "Test event bin");
    append(static_cast<pfunction>(&TestGSPIObservation::test_event_cube), "Test event cube");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set SPI optimizer test methods
 ***************************************************************************/
void TestGSPIOptimize::set(void)
{
    // Set test name
    name("SPI optimizers");

    // Append tests to test suite
    append(static_cast<pfunction>(&TestGSPIOptimize::test_binned_optimizer), "Test binned optimizer");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GSPIInstDir class
 ***************************************************************************/
void TestGSPIResponse::test_inst_dir(void)
{
    // Test constructors
    test_try("Test constructors");
    try {
        // Void constructor
        GSPIInstDir dir1;

        // Copy constructor
        GSPIInstDir dir2(dir1);

        // If we arrived here, signal success
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Create void object
    GSPIInstDir dir;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GSPIPointing class
 ***************************************************************************/
void TestGSPIResponse::test_pointing(void)
{
    // Test constructors
    test_try("Test constructors");
    try {
        // Void constructor
        GSPIPointing pnt1;

        // Copy constructor
        GSPIPointing pnt2(pnt1);

        // If we arrived here, signal success
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Create void object
    GSPIPointing pnt;

    // dir method
    GSkyDir sky;
    sky.radec_deg(37.0, 45.3);
    pnt.dir(sky);
    test_assert(pnt.dir() == sky, "Test dir() method.",
                "Expected "+sky.print()+", found "+pnt.dir().print());

    // copy constructor
    GSPIPointing pnt_copy(pnt);
    test_assert(pnt_copy.dir() == sky, "Test copy constructor method.",
                "Expected "+sky.print()+", found "+pnt_copy.dir().print());

    // assignment operator
    GSPIPointing pnt_assign = pnt;
    test_assert(pnt_assign.dir() == sky, "Test assignment operator method.",
                "Expected "+sky.print()+", found "+pnt_assign.dir().print());

    // clone
    GSPIPointing* pnt_clone = pnt.clone();
    test_assert(pnt_clone->dir() == sky, "Test clone() method.",
                "Expected "+sky.print()+", found "+pnt_clone->dir().print());

    // clear
    pnt.clear();
    sky.clear();
    test_assert(pnt.dir() == sky, "Test clear() method.",
                "Expected "+sky.print()+", found "+pnt.dir().print());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Checks handling of IAQ response files
 *
 * This function checks the handling of IAQ response files. IAQ response
 * files are 2D images that show the instrument response as function of
 * geometrical (Phi_geo) and measured (Phi_bar) Compton scatter angle.
 ***************************************************************************/
void TestGSPIResponse::test_response(void)
{
    // Test constructors
    test_try("Test constructors");
    try {
        // Void constructor
        GSPIResponse rsp1;

        // Copy constructor
        GSPIResponse rsp2(rsp1);

        // If we arrived here, signal success
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Checks handling of binned observations
 ***************************************************************************/
void TestGSPIObservation::test_binned_obs(void)
{
    // Test constructors
    test_try("Test constructors");
    try {
        // Void constructor
        GSPIObservation obs1;

        // Copy constructor
        GSPIObservation obs2(obs1);

        // If we arrived here, signal success
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Checks handling of SPI event bin
 ***************************************************************************/
void TestGSPIObservation::test_event_bin(void)
{
    // Test SPI event bin methods (one by one)
    test_try("Test event bin methods");
    try {
        // Event bin void constructor
        GSPIEventBin bin;

        // Copy constructor
        GSPIEventBin bin2(bin);

        // Assignment operator
        GSPIEventBin bin3 = bin;

        // clear method
        bin.clear();

        // clone method
        GSPIEventBin* bin4 = bin.clone();

        // size method
        test_value(bin.size(), 0.0, 1.0e-10, "Test size() method.");

        // dir method
        GSPIInstDir dir = bin.dir();

        // energy method
        GEnergy energy = bin.energy();

        // time method
        GTime time = bin.time();

        // counts method
        test_value(bin.counts(), 0.0, 1.0e-10, "Test counts() method for zero counts.");

        // error method
        test_value(bin.error(), 0.0, 1.0e-10, "Test error() method for zero counts.");

        // If we arrived here, signal success
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Checks handling of SPI event cube
 ***************************************************************************/
void TestGSPIObservation::test_event_cube(void)
{
    // Event cube void constructor
    GSPIEventCube cube;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test binned optimizer
 ***************************************************************************/
void TestGSPIOptimize::test_binned_optimizer(void)
{
    // Declare observations
    GObservations   obs;
    GSPIObservation run;

    // Exit test
    return;

}


/***************************************************************************
 * @brief Main entry point for test executable
 ***************************************************************************/
int main(void)
{
    // Allocate test suit container
    GTestSuites testsuites("SPI instrument specific class testing");

    // Set GAMMALIB_CALDB environment variable
    std::string caldb = "GAMMALIB_CALDB="+com_caldb;
    putenv((char*)caldb.c_str());

    // Initially assume that we pass all tests
    bool success = true;

    // Create test suites and append them to the container
    TestGSPIResponse    rsp;
    TestGSPIObservation obs;
    TestGSPIOptimize    opt;
    testsuites.append(rsp);
    testsuites.append(obs);
    testsuites.append(opt);

    // Run the testsuites
    success = testsuites.run();

    // Save test report
    testsuites.save("reports/GSPI.xml");

    // Return success status
    return (success ? 0 : 1);
}
