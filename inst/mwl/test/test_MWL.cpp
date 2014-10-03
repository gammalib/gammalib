/***************************************************************************
 *               test_MWL.cpp - Test multi-wavelength classes              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @brief Implementation of multi-wavelength classes unit tests
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <iostream>
#include "GMWLLib.hpp"
#include "test_MWL.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Constants __________________________________________________________ */
const std::string datadir        = PACKAGE_SOURCE"/inst/mwl/test/data";
const std::string lat_crab_model = datadir+"/crab.xml";
const std::string lat_crab_fits  = datadir+"/crab.fits";
const std::string crab_model     = datadir+"/crab_mwl.xml";
const std::string crab_fits      = datadir+"/crab_mwl.fits";
const std::string mwl_xml        = datadir+"/obs_mwl.xml";


/***********************************************************************//**
 * @brief Set test methods
 ***************************************************************************/
void TestGMWL::set(void)
{
    // Set test name
    name("GMWL");

    // Append tests to test suite
    append(static_cast<pfunction>(&TestGMWL::test_obs),       "Test observation handling");
    append(static_cast<pfunction>(&TestGMWL::test_optimizer), "Test optimizer");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGMWL* TestGMWL::clone(void) const
{
    // Clone test suite
    return new TestGMWL(*this);
}


/***********************************************************************//**
 * @brief Test observation handling
 ***************************************************************************/
void TestGMWL::test_obs(void)
{
    // Set filenames
    const std::string file1 = "test_mwl_obs.xml";

    // Construct observation
    test_try("Construct observation");
    try {
        GMWLObservation run1;
        GMWLObservation run2(lat_crab_fits);
        GMWLObservation run3 = run2;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Load observation
    test_try("Load observation");
    try {
        GMWLObservation run;
        run.load(lat_crab_fits);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Add observation to container
    test_try("Add observation to container");
    try {
        GMWLObservation run;
        run.load(lat_crab_fits);
        GObservations obs;
        run.id("001");
        obs.append(run);
        run.id("002");
        obs.append(run);
        run.id("003");
        obs.append(run);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test XML loading
    test_try("Test XML loading");
    try {
        GObservations obs = GObservations(mwl_xml);
        obs.save(file1);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Exit test
    return;
 
}


/***********************************************************************//**
 * @brief Test optimizer
 ***************************************************************************/
void TestGMWL::test_optimizer(void)
{
    // Declare observations
    GObservations   obs;

    // Load multi-wavelength observations
    test_try("Load multi-wavelength observations");
    try {
        GMWLObservation lat(lat_crab_fits);
        obs.append(lat);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Setup model
    test_try("Setup model");
    try {
        GModels models;
        models.load(lat_crab_model);
        obs.models(models);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Perform LM optimization
    test_try("Perform LM optimization");
    try {
        GLog log;
        log.cout(false);
        GOptimizerLM opt(log);
        opt.max_iter(1000);
        obs.optimize(opt);
        obs.errors(opt);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Load multi-wavelength observations
    test_try("Load multi-wavelength observations");
    try {
        obs.clear();
        GMWLObservation comptel(crab_fits, "COMPTEL");
        GMWLObservation lat(crab_fits, "LAT");
        GMWLObservation hess(crab_fits, "HESS");
        obs.append(comptel);
        obs.append(lat);
        obs.append(hess);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Setup model
    test_try("Setup model");
    try {
        GModels models;
        models.load(crab_model);
        obs.models(models);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Perform LM optimization
    double fit_results[] = {83.6331,             0,
                            22.0145,             0,
                            0.00231438040883405, 0.000191759607698647,
                            -3.1111154491587,    0.125598,
                            1,                   0,
                            1,                   0,
                            83.6331,             0,
                            22.0145,             0,
                            4.5050091460345e-11, 1.95885803229878e-12,
                            -2.13740092737976,   0.00578641,
                            1000,                0,
                            1,                   0};
    test_try("Perform LM optimization");
    try {
        GLog log;
        log.cout(false);
        GOptimizerLM opt(log);
        opt.max_iter(1000);
        obs.optimize(opt);
        obs.errors(opt);
        test_try_success();
        for (int i = 0, j = 0; i < obs.models().size(); ++i) {
            const GModel* model = obs.models()[i];
            for (int k = 0; k < model->size(); ++k) {
                GModelPar par = (*model)[k];
                std::string msg = "Verify optimization result for " + par.print();
                test_value(par.value(), fit_results[j++], 1.0e-4, msg);
                test_value(par.error(), fit_results[j++], 1.0e-3, msg);
            }
        }
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Exit test
    return;

}


/***************************************************************************
 * @brief Main entry point for test executable
 ***************************************************************************/
int main(void)
{
    // Allocate test suit container
    GTestSuites testsuites("Multi-wavelength instrument specific class testing");

    // Initially assume that we pass all tests
    bool success = true;

    // Create a test suite
    TestGMWL test;

    // Append test suite to the container
    testsuites.append(test);

    // Run the testsuites
    success = testsuites.run();

    // Save test report
    testsuites.save("reports/GMWL.xml");

    // Return success status
    return (success ? 0 : 1);
}
