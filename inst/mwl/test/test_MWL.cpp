/***************************************************************************
 *              test_MWL.cpp  -  test multi-wavelength classes             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file test_MWL.cpp
 * @brief Testing of MWL classes.
 * @author J. Knodlseder
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
const std::string lat_crab_fits = "../inst/mwl/test/data/crab.fits";


/***********************************************************************//**
 * @brief Test observation handling
 ***************************************************************************/
void test_obs(void)
{
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
/*
    // Number of observations in data
    int nobs = 1;

    // Declare observations
    GObservations   obs;
    GCTAObservation run;

    // Load binned CTA observation
    try {
        run.load_binned(cta_cntmap);
        run.response(cta_irf,cta_caldb);
        obs.append(run);
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to load CTA run."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Setup GModels for optimizing
    GModels models = crab_plaw();
    obs.models(models);

    // Perform LM optimization
    GOptimizerLM opt;
    try {
        opt.max_iter(1000);
        obs.optimize(opt);
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to perform LM optimization."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";
    std::cout << opt;
    std::cout << obs << std::endl;
*/
    // Plot final test success
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
