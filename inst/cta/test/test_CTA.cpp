/***************************************************************************
 *                      test_CTA.cpp  -  test CTA classes                  *
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
 * @file test_CTA.cpp
 * @brief Testing of CTA classes.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <iostream>
#include "GCTALib.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Constants __________________________________________________________ */
const std::string cta_caldb  = "../inst/cta/test/irf";
const std::string cta_irf    = "kb_E_50h_v3";
const std::string cta_events = "../inst/cta/test/data/run_00006028_eventlist_reco.fits.gz";
const double twopi           =  6.283185307179586476925286766559005768394;


/***********************************************************************//**
 * @brief Setup Crab point source powerlaw model.
 ***************************************************************************/
GModels crab_plaw(void)
{
    // Setup GModels for optimizing
    GModelSpatialPtsrc point_source;
    GModelSpectralPlaw power_law;
    GModel             crab;
    GModels            models;
    try {
        GSkyDir dir;
        //dir.radec_deg(83.6331, +22.0145);
        dir.radec_deg(117.0, -33.0);  // Adapt to source position in file
        point_source = GModelSpatialPtsrc(dir);
        power_law    = GModelSpectralPlaw(1.0e-7, -2.1);
        crab         = GModel(point_source, power_law);
        crab.name("Crab");
        models.append(crab);
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable setup GModels for optimizing." 
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    
    // Return model
    return models;
}


/***********************************************************************//**
 * @brief Test CTA response handling.
 ***************************************************************************/
void test_response(void)
{
    // Dump header
    std::cout << "Test GCTAResponse: ";

    // Test CTA response loading
    try {
        // Load response
        GCTAResponse rsp;
        rsp.caldb(cta_caldb);
        rsp.load(cta_irf);
        //std::cout << rsp << std::endl;
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to load CTA response." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test CTA Aeff response
    try {
        // Load response
        GCTAResponse rsp;
        rsp.caldb(cta_caldb);
        rsp.load(cta_irf);

        // Sum over effective area for control
        GEnergy      eng;
        GCTAInstDir  obsDir;
        GSkyDir      srcDir;
        GTime        time;
        GCTAPointing pnt;
        double sum = 0.0;
        double ref = 123299125000.0;
        for (int i = 0; i < 30; ++i) {
            eng.TeV(pow(10.0, -1.7 + 0.1*double(i)));
            double aeff = rsp.aeff(obsDir, eng, time, srcDir, eng, time, pnt);
            //std::cout << eng << " " << aeff << std::endl;
            sum += aeff;
        }
        //printf("%30.5f\n", sum);
        if (fabs(sum - ref) > 0.1) {
            std::cout << std::endl
                      << "TEST ERROR: Effective area differs from expected value"
                      << " (difference=" << fabs(sum - ref) << ")"
                      << std::endl;
            throw;
        }
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Error occured in CTA Aeff response." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test CTA Psf response
    try {
        // Load response
        GCTAResponse rsp;
        rsp.caldb(cta_caldb);
        rsp.load(cta_irf);
        
        // Integrate Psf
        GEnergy      eng;
        GTime        time;
        GCTAInstDir  obsDir;
        GSkyDir      srcDir;
        GCTAPointing pnt;
        srcDir.radec_deg(0.0, 0.0);
        for (double e = 0.1; e < 10.0; e*=2) {
            eng.TeV(e);
            double r     = 0.0;
            double dr    = 0.001;
            int    steps = 1.0/dr;
            double sum   = 0.0;
            for (int i = 0; i < steps; ++i) {
                r   += dr;
                obsDir.radec(0.0, r);
                sum += rsp.psf(obsDir, eng, time, srcDir, eng, time, pnt) * twopi * r * dr;
            }
            if ((sum - 1.0) > 0.001) {
                std::cout << std::endl
                        << "TEST ERROR: Psf integral differs from expected value"
                        << " (difference=" << (sum - 1.0) << ")"
                        << std::endl;
                throw;
            }
        }
        
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Error occured in CTA Psf response." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Dump final ok
    std::cout << " ok." << std::endl;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test unbinned observation handling.
 ***************************************************************************/
void test_unbinned_obs(void)
{
    // Write header
    std::cout << "Test CTA unbinned observation handling: ";

    // Declare observations
    GObservations   obs;
    GCTAObservation run;

    // Load unbinned CTA observation
    try {
        run.load_unbinned(cta_events);
        run.response(cta_irf,cta_caldb);
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to load CTA run."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Add observation (twice) to data
    try {
        obs.append(run);
        obs.append(run);
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to append CTA run to observations." 
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Loop over all events using iterators
    try {
        int num = 0;
        for (GObservations::iterator event = obs.begin(); event != obs.end(); ++event) {
            //std::cout << *event->energy() << std::endl;
            //std::cout << event->test() << std::endl;
            num++;
        }
        if (num != 17020) {
            std::cout << std::endl
                      << "TEST ERROR: Wrong number of iterations in GObservations::iterator."
                      << " (excepted 17020, found " << num << ")" << std::endl;
            throw;
        }
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to iterate GObservations." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Loop over all events using iterator
    try {
        int num = 0;
        GCTAEventList *ptr = (GCTAEventList*)run.events();
        for (GCTAEventList::iterator event = ptr->begin(); event != ptr->end(); ++event) {
            //std::cout << *event->energy() << std::endl;
            num++;
        }
        if (num != 8510) {
            std::cout << std::endl
                      << "TEST ERROR: Wrong number of iterations in GCTAEventList::iterator."
                      << " (excepted 8510, found " << num << ")" << std::endl;
            throw;
        }
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to iterate GCTAEventList." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Plot final test success
    std::cout << " ok." << std::endl;

    // Exit test
    return;
 
}


/***********************************************************************//**
 * @brief Test unbinned optimizer.
 ***************************************************************************/
void test_unbinned_optimizer(void)
{
    // Write header
    std::cout << "Test unbinned optimizer: ";
    
    // Number of observations in data
    int nobs = 1;

    // Declare observations
    GObservations   obs;
    GCTAObservation run;

    // Load unbinned CTA observation
    try {
        run.load_unbinned(cta_events);
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

    // Setup LM optimizer
    GOptimizerLM opt;
    try {
        obs.optimize(opt);
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable setup optimizer."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";
    std::cout << obs << std::endl;

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
    std::cout << "*****************************************" << std::endl;
    std::cout << "* CTA instrument specific class testing *" << std::endl;
    std::cout << "*****************************************" << std::endl;

    // Execute the tests
    test_response();
    test_unbinned_obs();
    test_unbinned_optimizer();

    // Return
    return 0;
}
