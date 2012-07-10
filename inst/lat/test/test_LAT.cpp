/***************************************************************************
 *                      test_LAT.cpp  -  test LAT classes                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
 * @file test_LAT.cpp
 * @brief Testing of LAT classes.
 * @author J. Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <iostream>
#include <unistd.h>
#include "GLATLib.hpp"
#include "GTools.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Constants __________________________________________________________ */
const std::string lat_caldb = "../inst/lat/caldb";
const std::string dirPass6  = "../inst/lat/test/data/p6v3";
const std::string dirPass7  = "../inst/lat/test/data/p7v6";


/***********************************************************************//**
 * @brief Livetime cube functions
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
 * @brief Test one specific response
 *
 * @param[in] irf Instrument response function.
 *
 * Verifies the ability to load and to save Fermi/LAT response functions.
 ***************************************************************************/
void test_one_response(const std::string& irf)
{
    // Set FITS filename
    std::string fitsfile = "test_rsp_" + irf + ".fits";
    
    // Remove FITS file
    std::string cmd = "rm -rf " + fitsfile;
    system(cmd.c_str());

    // Try loading the response
    try {
        GLATResponse rsp;
        rsp.caldb(lat_caldb);
        std::cout << ".";
        rsp.load(irf+"::front");
        std::cout << ".";
        rsp.load(irf+"::back");
        std::cout << ".";
        rsp.load(irf);
        //std::cout << std::endl << rsp << std::endl;
        std::cout << ".";
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to load LAT response "
                  << irf << "." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Try saving the response
    try {
        GLATResponse rsp;
        rsp.caldb(lat_caldb);
        rsp.load(irf);
        rsp.save(fitsfile);
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to save LAT response "
                  << irf << "." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test livetime cube handling for a specific dataset
 *
 * @param[in] datadir Directory of test data.
 *
 * Verifies the ability to handle Fermi/LAT livetime cubes.
 ***************************************************************************/
void test_one_ltcube(const std::string& datadir, const double& reference)
{
    // Set filenames
    std::string lat_ltcube = datadir+"/ltcube.fits";
    std::string file1      = "test_lat_ltcube.fits";
    std::string file2      = "test_lat_ltcube_phi.fits";

    // Load livetime cube
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

    // Test operators (no phi dependence)
    try {
        // Load livetime cube
        GLATLtCube ltcube(lat_ltcube);

        // Initialise sky direction and energy
        GSkyDir dir;
        GEnergy energy;

        // Test cos theta integration operator.
        double sum = ltcube(dir, energy, test_fct1);
        if (fabs(sum-reference) > 0.001) {
            std::cout << std::endl
                      << "TEST ERROR: Livetime cube sum "
                      << sum << " is invalid, "
                      << reference << " expected."
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
        GLATLtCube ltcube(lat_ltcube);

        // Initialise sky direction and energy
        GSkyDir dir;
        GEnergy energy;

        // Test cos theta and phi integration operator. The sum differs from
        // above since the actual test dataset does not cover the same time
        // interval
        double sum = ltcube(dir, energy, test_fct2);
        if (fabs(sum-0.0) > 0.001) {
            std::cout << std::endl
                      << "TEST ERROR: Livetime cube sum "
                      << sum << " is invalid, "
                      << 0.0 << " expected."
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
        GLATLtCube ltcube(lat_ltcube);

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
        GLATLtCube ltcube(lat_ltcube);

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

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test unbinned observation handling for a specific dataset
 *
 * @param[in] datadir Directory of test data.
 *
 * Verifies the ability to handle unbinned Fermi/LAT data.
 ***************************************************************************/
void test_one_unbinned_obs(const std::string& datadir)
{
    // Set filenames
    std::string lat_ft1       = datadir+"/ft1.fits";
    std::string lat_ft2       = datadir+"/ft2.fits";
    std::string lat_unbin_xml = datadir+"/obs_unbinned.xml";
    std::string file1         = "test_lat_obs_unbinned.xml";

    // Declare observations
    GObservations   obs;
    GLATObservation run;

    // Determine number of events in FT1 file
    GFits ft1(lat_ft1);
    int nevents = ft1.table("EVENTS")->nrows();
    ft1.close();

    // Load unbinned LAT observation
    try {
        run.load_unbinned(lat_ft1, lat_ft2, "");
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to load unbinned LAT data"
                  << " (" << datadir << ")."
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
                  << "TEST ERROR: Unable to add LAT data to observations"
                  << " (" << datadir << ")." 
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
        if (num != nevents*2) {
            std::cout << std::endl
                      << "TEST ERROR: " << num
                      << " iterations in GObservations::iterator instead of"
                      << " expected " << nevents*2 << " iterations."
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
        if (num != nevents) {
            std::cout << std::endl
                      << "TEST ERROR: " << num
                      << " iterations in GLATEventList::iterator instead of "
                      << " expected " << nevents << " iterations."
                      << std::endl;
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

    // Test XML loading
    try {
        // Load observations
        obs = GObservations(lat_unbin_xml);
        
        // Save observations
        obs.save(file1);
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to load unbinned observation from XML file."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Exit test
    return;

}


/***********************************************************************//**
 * @brief Test binned observation handling for a specific dataset
 *
 * @param[in] datadir Directory of test data.
 *
 * Verifies the ability to handle binned Fermi/LAT data.
 ***************************************************************************/
void test_one_binned_obs(const std::string& datadir, const std::string& irf)
{
    // Set filenames
    std::string lat_cntmap  = datadir+"/cntmap.fits";
    std::string lat_srcmap  = datadir+"/srcmap.fits";
    std::string lat_expmap  = datadir+"/binned_expmap.fits";
    std::string lat_ltcube  = datadir+"/ltcube.fits";
    std::string lat_bin_xml = datadir+"/obs_binned.xml";
    std::string file1       = "test_lat_obs_binned.xml";

    // Declare observations
    GObservations   obs;
    GLATObservation run;

    // Determine number of bins and events in counts map
    GFits       cntmap(lat_cntmap);
    GFitsImage* image   = cntmap.image(0);
    double      nevents = 0.0;
    int         nsize   = image->size();
    for (int i = 0; i < nsize; ++i) {
        nevents += image->pixel(i);
    }
    cntmap.close();

    // Load LAT binned observation from counts map
    try {
        run.load_binned(lat_cntmap, "", "");
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to load binned LAT data from cntmap"
                  << " (" << datadir << ")."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Reload LAT binned observation from source map
    try {
        run.load_binned(lat_srcmap, "", "");
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to load binned LAT data from srcmap"
                  << " (" << datadir << ")."
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
                  << "TEST ERROR: Unable to add LAT data to observations"
                  << " (" << datadir << ")." 
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
        if (sum != 2*nevents) {
            std::cout << std::endl
                      << "TEST ERROR: " << sum
                      << " events in GObservations instead of"
                      << " expected " << 2*nevents << " events."
                      << std::endl;
            throw;
        }
        if (num != 2*nsize) {
            std::cout << std::endl
                      << "TEST ERROR: " << num
                      << " iterations in GObservations::iterator instead of "
                      << " expected " << 2*nsize << " iterations."
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
            num++;
            sum += (int)event->counts();
        }
        if (sum != nevents) {
            std::cout << std::endl
                      << "TEST ERROR: " << sum
                      << " events in GLATEventCube instead of"
                      << " expected " << 2*nevents << " events."
                      << std::endl;
            throw;
        }
        if (num != nsize) {
            std::cout << std::endl
                      << "TEST ERROR: " << num
                      << " iterations in GLATEventCube::iterator instead of "
                      << " expected " << 2*nsize << " iterations."
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
        run.response(irf, lat_caldb);

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

    // Test XML loading
    try {
        // Load observations
        obs = GObservations(lat_bin_xml);
        
        // Save observations
        obs.save(file1);
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to load binned observation from XML file."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Exit test
    return;

}


/***********************************************************************//**
 * @brief Test binned optimizer
 *
 * @param[in] datadir Directory of test data.
 *
 * Verifies the ability optimize binned Fermi/LAT data.
 ***************************************************************************/
void test_one_binned_optimizer(const std::string& datadir, const std::string& irf)
{
    // Set filenames
    std::string lat_srcmap    = datadir+"/srcmap.fits";
    std::string lat_expmap    = datadir+"/binned_expmap.fits";
    std::string lat_ltcube    = datadir+"/ltcube.fits";
    std::string lat_model_xml = datadir+"/source_model.xml";

    // Setup GObservations for optimizing
    GObservations   obs;
    GLATObservation run;
    try {
        run.load_binned(lat_srcmap, lat_expmap, lat_ltcube);
        run.response(irf, lat_caldb);
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
    obs.models(lat_model_xml);

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
    std::cout << obs.models() << std::endl;

    // Exit test
    return;

}


/***********************************************************************//**
 * @brief Test Fermi/LAT response handling
 *
 * Verifies all Fermi/LAT responses.
 ***************************************************************************/
void test_response(void)
{
    // Write header
    std::cout << "Test response: ";

    // Test various IRFs
    test_one_response("P6_v3_diff");
    test_one_response("P7SOURCE_V6");

    // Plot final test success
    std::cout << " ok." << std::endl;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test livetime cube handling
 *
 * Verifies handling of Fermi/LAT livetime cubes.
 ***************************************************************************/
void test_ltcube(void)
{
    // Write header
    std::cout << "Test livetime cube: ";

    // Test various datasets
    test_one_ltcube(dirPass6, 339733.528629);
    test_one_ltcube(dirPass7, 248009.734604);

    // Plot final test success
    std::cout << " ok." << std::endl;

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test unbinned observation handling
 *
 * Verifies handling of unbinned data. 
 ***************************************************************************/
void test_unbinned_obs(void)
{
    // Write header
    std::cout << "Test unbinned observation handling: ";

    // Test various datasets
    test_one_unbinned_obs(dirPass6);
    test_one_unbinned_obs(dirPass7);

    // Plot final test success
    std::cout << " ok." << std::endl;

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test binned observation handling
 *
 * Verifies the ability to handle binned Fermi/LAT data.
 ***************************************************************************/
void test_binned_obs(void)
{
    // Write header
    std::cout << "Test binned observation handling: ";

    // Test various datasets
    test_one_binned_obs(dirPass6, "P6_v3_diff");
    test_one_binned_obs(dirPass7, "P7SOURCE_V6");

    // Plot final test success
    std::cout << " ok." << std::endl;

    // Exit test
    return;

}


/***********************************************************************//**
 * @brief Test binned optimizer handling
 *
 * Verifies the ability to handle binned Fermi/LAT optimization.
 ***************************************************************************/
void test_binned_optimizer(void)
{
    // Write header
    std::cout << "Test binned optimizer: ";

    // Test various datasets
    test_one_binned_optimizer(dirPass6, "P6_v3_diff");
    test_one_binned_optimizer(dirPass7, "P7SOURCE_V6");

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

    // Check if data directory exists
    bool has_data = (access("../inst/lat/test/data", R_OK) == 0);

    // Execute tests that are not requiring data
    test_response();

    // Execute tests requiring data
    if (has_data) {

        // Set CALDB environment variable
        std::string caldb = "CALDB="+lat_caldb;
        putenv((char*)caldb.c_str());

        // Execute tests
        test_ltcube();
        test_unbinned_obs();
        test_binned_obs();
        test_binned_optimizer();
    }
    else {
        std::cout << "Skipped several tests since no test data have been found." << std::endl;
    }

    // Return
    return 0;
}
