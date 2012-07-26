/***************************************************************************
 *                      test_COM.cpp  -  test COM classes                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by <author>                                         *
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
 * @file test_COM.cpp
 * @brief Testing of COM classes
 * @author <author>
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <iostream>
#include <unistd.h>
#include "GCOMLib.hpp"
#include "GTools.hpp"

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
const std::string com_obs       = datadir+"obs.xml";


/***********************************************************************//**
 * @brief Checks handling of IAQ response files
 *
 * This function checks the handling of IAQ response files. IAQ response
 * files are 2D images that show the instrument response as function of
 * geometrical (Phi_geo) and measured (Phi_bar) Compton scatter angle.
 ***************************************************************************/
void test_iaq_response(void)
{
    // Test response loading
    try {
        // Construct observation from datasets
        GCOMResponse rsp(com_caldb, com_iaq);
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to construct IAQ response."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Checks handling of COMPTEL response files
 ***************************************************************************/
void test_response(void)
{
    // Dump header
    std::cout << "Test COMPTEL response: ";

    // Test IAQ
    test_iaq_response();

    // Dump final ok
    std::cout << " ok." << std::endl;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Checks handling of binned observations
 *
 * This function checks the handling of binned observations. Binned
 * observations are defined by
 * a DRE file containing the events,
 * a DRB file containing the background model,
 * a DRG file containing the geometrical probability of having an interaction
 *            in the second detector plane, and
 * a DRX file containing the sky exposure.
 *
 * Note that in first approximation we can use the DRG file as background
 * model.
 ***************************************************************************/
void test_binned_obs(void)
{
    // Write header
    std::cout << "Test binned COMPTEL observation handling: ";

    // Declare observation and container
    GObservations obs;

    // Dataset constructor
    try {
        // Construct observation from datasets
        GCOMObservation com(com_dre, com_drb, com_drg, com_drx);
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to construct observation from datasets."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // XML constructor
    try {
        // Construct observation from XML file
        GObservations obs(com_obs);
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to construct observation from XML file."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Notify final test success
    std::cout << " ok." << std::endl;

    // Exit test
    return;
 
}


/***********************************************************************//**
 * @brief Main test function
 ***************************************************************************/
int main(void)
{
    // Dump header
    std::cout << std::endl;
    std::cout << "*********************************************" << std::endl;
    std::cout << "* COMPTEL instrument specific class testing *" << std::endl;
    std::cout << "*********************************************" << std::endl;

    // Check if data directory exists
    bool has_data = (access(datadir.c_str(), R_OK) == 0);

    // Execute tests not needing data
    test_response();
    
    // Execute tests requiring data
    if (has_data) {

        // Set CALDB environment variable
        std::string caldb = "CALDB="+com_caldb;
        putenv((char*)caldb.c_str());

        // Execute tests
        test_binned_obs();
    }
    else {
        std::cout << "Skipped several tests since no test data have been found."
                  << std::endl;
    }

    // Return
    return 0;
}
