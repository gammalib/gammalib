/***************************************************************************
 *                      test_LAT.cpp  -  test LAT classes                  *
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
 * @file test_LAT.cpp
 * @brief Testing of LAT classes.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <iostream>
#include "GLATLib.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Constants __________________________________________________________ */
const std::string lat_caldb  = "test/irf";
const std::string lat_irf    = "Pass5_v0";
const std::string lat_ft1    = "test/data/ft1.fits.gz";
const std::string lat_ft2    = "test/data/ft2.fits.gz";
const std::string lat_cntmap = "test/data/cntmap.fits.gz";
const std::string lat_srcmap = "test/data/srcmap.fits.gz";
const double      twopi      =  6.283185307179586476925286766559005768394;


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
        dir.radec_deg(83.6331, +22.0145);
        point_source = GModelSpatialPtsrc(dir);
        power_law    = GModelSpectralPlaw(1.0e-7, -2.1);
        crab         = GModel(point_source, power_law);
        crab.name("Crab");
        models.add(crab);
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
 * @brief Test response handling.
 ***************************************************************************/
void test_response(void)
{
    // Remove FITS file
    system("rm -rf test_rsp.fits");

    std::cout << "Test GLATResponse: ";
    try {
        // Load response
        GLATResponse rsp;
        rsp.caldb(lat_caldb);
        rsp.load(lat_irf, "front");

        // Save response
        rsp.save("test_rsp.fits");
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to load LAT response." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ". ok." << std::endl;
}


/***********************************************************************//**
 * @brief Test unbinned observation handling.
 ***************************************************************************/
void test_unbinned_obs(void)
{
    // Write header
    std::cout << "Test unbinned observation handling: ";

    // Declare observations
    GObservations   data;
    GLATObservation obs;

    // Load unbinned LAT observation
    try {
        obs.load_unbinned(lat_ft1, lat_ft2, "");
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to load LAT observation."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Add observation (twice) to data
    try {
        data.append(obs);
        data.append(obs);
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to add LAT observation to data." 
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Loop over all events in GObservations using iterators
    try {
        int num = 0;
        for (GObservations::iterator event = data.begin(); event != data.end(); ++event) {
//            std::cout << *event << std::endl;
//            std::cout << event->test() << std::endl;
            num++;
        }
        if (num != 4038) {
            std::cout << std::endl << 
                      "TEST ERROR: Wrong number of iterations in GObservations::iterator."
                      << std::endl;
            throw;
        }
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to iterate GObservations." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Loop over all events using iterator
    try {
        int num = 0;
        GLATEventList *ptr = (GLATEventList*)obs.events();
        for (GLATEventList::iterator event = ptr->begin(); event != ptr->end(); ++event) {
//            std::cout << *event << std::endl;
            num++;
        }
        if (num != 2019) {
            std::cout << std::endl << 
                      "TEST ERROR: Wrong number of iterations in GLATEventList::iterator."
                      << std::endl;
            throw;
        }
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to iterate GLATEventList." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Plot final test success
    std::cout << " ok." << std::endl;

//GLATEventList *ptr = (GLATEventList*)obs.events();
//std::cout << obs << std::endl;
//std::cout << *ptr << std::endl;

    // Exit test
    return;
 
}


/***********************************************************************//**
 * @brief Test binned observation handling.
 ***************************************************************************/
void test_binned_obs(void)
{
    // Write header
    std::cout << "Test LAT binned observation handling: ";

    // Declare observations
    GObservations   data;
    GLATObservation obs;

    // Load LAT binned observation
    try {
        obs.load_binned(lat_cntmap, "", "");
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to load binned LAT observation." 
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Reload LAT binned observation from srcmap
    try {
        obs.load_binned(lat_srcmap, "", "");
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to reload binned LAT observation." 
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Add observation (twice) to data
    try {
        data.append(obs);
        data.append(obs);
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to add LAT observation to data." 
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Loop over all events in GObservations using iterators
    try {
        int num = 0;
        int sum = 0;
        for (GObservations::iterator event = data.begin(); event != data.end(); ++event) {
            num++;
            sum += (int)event->counts();
        }
        if (sum != 2718 || num != 400000) {
            std::cout << std::endl << 
                      "TEST ERROR: Wrong number of iterations in GObservations::iterator."
                      << std::endl;
            throw;
        }
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to iterate GObservations." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Loop over all events using iterator
    try {
        int num = 0;
        int sum = 0;
        GLATEventCube *ptr = (GLATEventCube*)obs.events();
        for (GLATEventCube::iterator event = ptr->begin(); event != ptr->end(); ++event) {
//            std::cout << *((GLATEventBin*)&(*event));
            num++;
            sum += (int)event->counts();
        }
        if (sum != 1359 || num != 200000) {
            std::cout << std::endl << 
                      "TEST ERROR: Wrong number of iterations in GLATEventCube::iterator."
                      << std::endl;
            throw;
        }
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to iterate GLATEventCube." << std::endl;
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
 * @brief Test binned optimizer.
 ***************************************************************************/
void test_binned_optimizer(void)
{
    // Write header
    std::cout << "Test binned optimizer: ";
    
    // Number of observations in data
    int nobs = 1;

    // Setup GObservations for optimizing
    GObservations   data;
    GLATObservation obs;
    try {
        obs.load_binned(lat_cntmap, "", "");
        for (int i = 0; i < nobs; ++i)
            data.append(obs);
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable setup GData for optimizing." 
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";
//std::cout << data << std::endl;

    // Setup GModels for optimizing
    GModels models = crab_plaw();
    data.models(models);

    // Setup parameters for optimizing
    GOptimizerPars pars;
    try {
        pars = GOptimizerPars(models);
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable setup GOptimizerPars for optimizing." 
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";
//std::cout << pars << std::endl;

    // Time simple GObservations iterator
    double  t_elapse;
    try {
        clock_t t_start = clock();
        int num = 0;
        int sum = 0;
        for (GObservations::iterator event = data.begin(); event != data.end(); ++event) {
            num++;
            sum += (int)event->counts();
        }
        t_elapse = double(clock() - t_start)/double(CLOCKS_PER_SEC);
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to iterate GObservations." 
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << "." << std::endl;
    std::cout << " - Reference time for GData::iterator: " << t_elapse << std::endl;

    // Setup LM optimizer
    GOptimizerLM opt;
    try {
        data.optimize(opt);
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable setup optimizer." 
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";
//std::cout << opt << std::endl;
std::cout << data << std::endl;

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
    std::cout << "* LAT instrument specific class testing *" << std::endl;
    std::cout << "*****************************************" << std::endl;

    // Execute the tests
    test_response();
    test_unbinned_obs();
    test_binned_obs();
    test_binned_optimizer();

    // Return
    return 0;
}
