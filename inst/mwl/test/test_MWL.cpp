/***************************************************************************
 *              test_MWL.cpp  -  test multi-wavelength classes             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Juergen Knoedlseder                         *
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
 * @file test_MWL.cpp
 * @brief Testing of MWL classes.
 * @author J. Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <iostream>
#include "GMWLLib.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Constants __________________________________________________________ */
const std::string datadir        = "../inst/mwl/test/data";
const std::string lat_crab_model = datadir+"/crab.xml";
const std::string lat_crab_fits  = datadir+"/crab.fits";
const std::string crab_model     = datadir+"/crab_mwl.xml";
const std::string crab_fits      = datadir+"/crab_mwl.fits";
const std::string mwl_xml        = datadir+"/obs_mwl.xml";


/***********************************************************************//**
 * @brief Test observation handling
 ***************************************************************************/
void test_obs(void)
{
    // Set filenames
    const std::string file1 = "test_mwl_obs.xml";

    // Write header
    std::cout << "Test observation handling: ";

    // Construct observation
    try {
        GMWLObservation run1;
        GMWLObservation run2(lat_crab_fits);
        GMWLObservation run3 = run2;
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to construct GMWLObservation."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Load observation
    try {
        GMWLObservation run;
        run.load(lat_crab_fits);
        //std::cout << run << std::endl;
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to load GMWLObservation."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Add observation to container
    try {
        GMWLObservation run;
        run.load(lat_crab_fits);
        GObservations obs;
        obs.append(run);
        obs.append(run);
        obs.append(run);
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to append GMWLObservation to GObservations."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test XML loading
    try {
        // Load observations
        GObservations obs = GObservations(mwl_xml);
        
        // Save observations
        obs.save(file1);
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to load MWL observation from XML file."
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
 * @brief Test binned optimizer.
 ***************************************************************************/
void test_optimizer(void)
{
    // Write header
    std::cout << "Test optimizer: ";

    // Declare observations
    GObservations   obs;

    // Load multi-wavelength observations
    try {
        GMWLObservation lat(lat_crab_fits);
        obs.append(lat);
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to append MWL observation(s) to container."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Setup model
    try {
        GModels models;
        models.load(lat_crab_model);
        obs.models(models);
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to setup model from XML file for fitting."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Perform LM optimization
    try {
        GLog log;
        log.cout(false);
        GOptimizerLM opt(log);
        opt.max_iter(1000);
        obs.optimize(opt);
        //std::cout << obs << std::endl;
        //std::cout << std::endl << opt << std::endl;
        //std::cout << *(obs.models()) << std::endl;
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to perform LM optimization."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Load multi-wavelength observations
    try {
        obs.clear();
        GMWLObservation comptel(crab_fits, "COMPTEL");
        GMWLObservation lat(crab_fits, "LAT");
        GMWLObservation hess(crab_fits, "HESS");
        obs.append(comptel);
        obs.append(lat);
        obs.append(hess);
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to append MWL observation(s) to container."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Setup model
    try {
        GModels models;
        models.load(crab_model);
        obs.models(models);
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to setup model from XML file for fitting."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Perform LM optimization
    try {
        GLog log;
        log.cout(false);
        GOptimizerLM opt(log);
        opt.max_iter(1000);
        obs.optimize(opt);
        std::cout << std::endl << obs << std::endl;
        std::cout << opt << std::endl;
        std::cout << obs.models() << std::endl;
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to perform LM optimization."
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
 * @brief Main test function .
 ***************************************************************************/
int main(void)
{
    // Dump header
    std::cout << std::endl;
    std::cout << "**********************************" << std::endl;
    std::cout << "* Multi-wavelength class testing *" << std::endl;
    std::cout << "**********************************" << std::endl;

    // Execute the tests
    test_obs();
    test_optimizer();

    // Return
    return 0;
}
