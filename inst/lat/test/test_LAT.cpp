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
const std::string lat_caldb      = "../inst/lat/test/irf";
//const std::string lat_irf        = "Pass5_v0";
const std::string lat_irf        = "P6_v3_diff";
const std::string lat_ft1        = "../inst/lat/test/data/ft1.fits";
const std::string lat_ft2        = "../inst/lat/test/data/ft2.fits";
const std::string lat_cntmap     = "../inst/lat/test/data/cntmap.fits";
const std::string lat_srcmap     = "../inst/lat/test/data/srcmap.fits";
const std::string lat_expmap     = "../inst/lat/test/data/binned_expmap.fits";
const std::string lat_ltcube     = "../inst/lat/test/data/ltcube.fits";
const std::string lat_ltcube_phi = "../inst/lat/test/data/ltcube_phi02.fits";
const std::string lat_xml        = "../inst/lat/test/data/source_model3.xml";
const double      twopi          =  6.283185307179586476925286766559005768394;


/***********************************************************************//**
 * @brief Test response handling.
 ***************************************************************************/
void test_response(void)
{
    // Remove FITS file
    system("rm -rf test_rsp.fits");

    // Write header
    std::cout << "Test response: ";

    // Try loading
    try {
        GLATResponse rsp;
        rsp.caldb(lat_caldb);
        std::cout << ".";
        rsp.load(lat_irf+"::front");
        std::cout << ".";
        rsp.load(lat_irf+"::back");
        std::cout << ".";
        rsp.load(lat_irf);
std::cout << std::endl << rsp << std::endl;
        std::cout << ".";
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to load LAT response." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Try saving
    try {
        GLATResponse rsp;
        rsp.caldb(lat_caldb);
        rsp.load(lat_irf);
        rsp.save("test_rsp.fits");
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to save LAT response." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Plot final test success
    std::cout << " ok." << std::endl;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Livetime cube functions.
 ***************************************************************************/
double test_fct1(const double& ctheta)
{
    //std::cout << ctheta << " ";
    return 1.0;
}
double test_fct2(const double& ctheta, const double& phi)
{
    //std::cout << "(" << ctheta << "," << phi << ")" << " ";
    return 1.0;
}


/***********************************************************************//**
 * @brief Test livetime cube handling.
 ***************************************************************************/
void test_ltcube(void)
{
    // Set filenames
    const std::string file1 = "test_lat_ltcube.fits";
    const std::string file2 = "test_lat_ltcube_phi.fits";

    // Write header
    std::cout << "Test livetime cube: ";

    // Load livetime cube (no phi dependence)
    try {
        // Load livetime cube
        GLATLtCube ltcube(lat_ltcube);
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to load LAT livetime cube."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Load livetime cube (no phi dependence)
    try {
        // Load livetime cube
        GLATLtCube ltcube(lat_ltcube_phi);
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to load phi-dependent LAT livetime cube."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test operators (no phi dependence)
    try {
        // Load livetime cube
        //GLATLtCube ltcube(lat_ltcube);
        GLATLtCube ltcube(lat_ltcube_phi);

        // Initialise sky direction and energy
        GSkyDir dir;
        GEnergy energy;

        // Test cos theta integration operator. The total ontime of the dataset is
        // 811432.0 sec, hence the efficiency is 42%. This seems reasonable as the
        // test function returns 1 for half of the sky.
        double sum = ltcube(dir, energy, test_fct1);
//        if (fabs(sum-247712.0724) > 0.001) {
        if (fabs(sum-339334.6556) > 0.001) {
            std::cout << std::endl
                      << "TEST ERROR: Invalid livetime cube sum (expected 339334.6556,"
                      << " encountered difference " << sum-339334.6556
                      << std::endl;
            throw;
        }
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to access LAT livetime cube."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test operators (phi dependence)
    try {
        // Load livetime cube
        GLATLtCube ltcube(lat_ltcube_phi);
        //std::cout << ltcube << std::endl;

        // Initialise sky direction and energy
        GSkyDir dir;
        GEnergy energy;

        // Test cos theta and phi integration operator. The sum differs from above
        // since the actual test dataset does not cover the same time interval
        double sum = ltcube(dir, energy, test_fct2);
        if (fabs(sum-339334.6641) > 0.001) {
            std::cout << std::endl
                      << "TEST ERROR: Invalid livetime cube sum (expected 339334.6641,"
                      << " encountered difference " << sum-339334.6641
                      << std::endl;
            throw;
        }
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to access phi-dependent LAT livetime cube."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Create livetime skymap (no phi dependence)
    try {
        // Allocate skymap
        GSkymap map("HPX", "GAL", 64, "RING", 1);
        
        // Load livetime cube
//        GLATLtCube ltcube(lat_ltcube);
        GLATLtCube ltcube(lat_ltcube_phi);

        // Create livetime skymap
        GEnergy energy;
        for (int i = 0; i < map.npix(); ++i) {
            GSkyDir dir = map.pix2dir(i);
            map(i) = ltcube(dir, energy, test_fct1);
        }

        // Save skymap
        map.save(file1, true);

    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to build LAT livetime cube skymap."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Create livetime skymap (phi dependence)
    try {
        // Allocate skymap
        GSkymap map("HPX", "GAL", 64, "RING", 1);
        
        // Load livetime cube
        GLATLtCube ltcube(lat_ltcube_phi);

        // Create livetime skymap
        GEnergy energy;
        for (int i = 0; i < map.npix(); ++i) {
            GSkyDir dir = map.pix2dir(i);
            map(i) = ltcube(dir, energy, test_fct2);
        }

        // Save skymap
        map.save(file2, true);

    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to build phi-dependent LAT livetime cube skymap."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Plot final test success
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
    std::cout << "Test unbinned observation handling: ";

    // Declare observations
    GObservations   obs;
    GLATObservation run;

    // Load unbinned LAT observation
    try {
        run.load_unbinned(lat_ft1, lat_ft2, "");
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to load unbinned LAT run."
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
        std::cout << std::endl << "TEST ERROR: Unable to add LAT run to observations." 
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Loop over all events in GObservations using iterators
    try {
        int num = 0;
        for (GObservations::iterator event = obs.begin(); event != obs.end(); ++event) {
            //std::cout << *(event->energy()) << std::endl;
            //std::cout << num << std::endl;
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
        std::cout << std::endl << "TEST ERROR: Unable to iterate GObservations."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Loop over all events using iterator
    try {
        int num = 0;
        GLATEventList *ptr = (GLATEventList*)run.events();
        for (GLATEventList::iterator event = ptr->begin(); event != ptr->end(); ++event) {
            //std::cout << *event->energy() << std::endl;
            num++;
        }
        if (num != 2019) {
            std::cout << std::endl <<
                      "TEST ERROR: Wrong number of iterations in GLATEventList::iterator."
                      << " (excepted 2019, found " << num << ")" << std::endl;
            throw;
        }
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to iterate GLATEventList."
                  << std::endl;
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
 * @brief Test binned observation handling.
 ***************************************************************************/
void test_binned_obs(void)
{
    // Write header
    std::cout << "Test binned observation handling: ";

    // Declare observations
    GObservations   obs;
    GLATObservation run;

    // Load LAT binned observation
    try {
        run.load_binned(lat_cntmap, "", "");
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to load binned LAT run."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Reload LAT binned observation from srcmap
    try {
        run.load_binned(lat_srcmap, "", "");
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to reload binned LAT run."
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
        std::cout << std::endl << "TEST ERROR: Unable to add LAT run to observation."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Loop over all events in GObservations using iterators
    try {
        int num = 0;
        int sum = 0;
        for (GObservations::iterator event = obs.begin(); event != obs.end(); ++event) {
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
        std::cout << std::endl << "TEST ERROR: Unable to iterate GObservations."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Loop over all events using iterator
    try {
        int num = 0;
        int sum = 0;
        GLATEventCube *ptr = (GLATEventCube*)run.events();
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
        std::cout << std::endl << "TEST ERROR: Unable to iterate GLATEventCube."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test mean PSF
    try {
        // Load obervation
        run.load_binned(lat_srcmap, lat_expmap, lat_ltcube);
        run.response(lat_irf,lat_caldb);

        // Initialise sky direction
        GSkyDir dir;

        // Set mean PSF
        GLATMeanPsf psf(dir, run);
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to setup mean LAT PSF."
                  << std::endl;
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

    // Setup GObservations for optimizing
    GObservations   obs;
    GLATObservation run;
    try {
        run.load_binned(lat_srcmap, lat_expmap, lat_ltcube);
        obs.append(run);
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to setup observation for optimizing."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Load models from XML file
    obs.models(lat_xml);

    // Setup LM optimizer
    GOptimizerLM opt;
    try {
        opt.max_iter(1000);
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
    std::cout << std::endl << opt << std::endl;
    std::cout << *obs.models() << std::endl;

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
    test_ltcube();
    test_unbinned_obs();
    test_binned_obs();
    test_binned_optimizer();

    // Return
    return 0;
}
