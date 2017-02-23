/***************************************************************************
 *                  test_LAT.cpp - test Fermi/LAT classes                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2016 by Juergen Knoedlseder                         *
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
 * @brief Testing of LAT classes
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <cstdlib>     // getenv
#include "GLATLib.hpp"
#include "GTools.hpp"
#include "test_LAT.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Constants __________________________________________________________ */
const std::string datadir   = std::getenv("TEST_LAT_DATA");
const std::string lat_caldb = datadir + "/../../caldb";
const std::string dirPass6  = datadir + "/p6v3";
const std::string dirPass7  = datadir + "/p7v6";
const std::string dirPass8  = datadir + "/p8v2";


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
 * @brief Set LAT response test methods
 ***************************************************************************/
void TestGLATResponse::set(void)
{
    // Set test name
    name("GLATResponse");

    // Append tests to test suite
    append(static_cast<pfunction>(&TestGLATResponse::test_response_p6),
           "Test P6 response");
    append(static_cast<pfunction>(&TestGLATResponse::test_response_p7),
           "Test P7 response");
    append(static_cast<pfunction>(&TestGLATResponse::test_response_p8),
           "Test P8 response");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGLATResponse* TestGLATResponse::clone(void) const
{
    // Clone test suite
    return new TestGLATResponse(*this);
}


/***********************************************************************//**
 * @brief Set LAT livetime cube test methods
 ***************************************************************************/
void TestGLATLtCube::set(void)
{
    // Set test name
    name("GLATLtCube");

    // Append tests to test suite
    append(static_cast<pfunction>(&TestGLATLtCube::test_ltcube_p6),
           "Test P6 livetime cube");
    append(static_cast<pfunction>(&TestGLATLtCube::test_ltcube_p7),
           "Test P7 livetime cube");
    append(static_cast<pfunction>(&TestGLATLtCube::test_ltcube_p8),
           "Test P8 livetime cube");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGLATLtCube* TestGLATLtCube::clone(void) const
{
    // Clone test suite
    return new TestGLATLtCube(*this);
}


/***********************************************************************//**
 * @brief Set LAT observation test methods
 ***************************************************************************/
void TestGLATObservation::set(void)
{
    // Set test name
    name("GLATObservation");

    // Append tests to test suite
    append(static_cast<pfunction>(&TestGLATObservation::test_unbinned_obs_p6),
           "Test P6 unbinned observation");
    append(static_cast<pfunction>(&TestGLATObservation::test_unbinned_obs_p7),
           "Test P7 unbinned observation");
    append(static_cast<pfunction>(&TestGLATObservation::test_unbinned_obs_p8),
           "Test P8 unbinned observation");
    append(static_cast<pfunction>(&TestGLATObservation::test_binned_obs_p6),
           "Test P6 binned observation");
    append(static_cast<pfunction>(&TestGLATObservation::test_binned_obs_p7),
           "Test P7 binned observation");
    append(static_cast<pfunction>(&TestGLATObservation::test_binned_obs_p8),
           "Test P8 binned observation");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGLATObservation* TestGLATObservation::clone(void) const
{
    // Clone test suite
    return new TestGLATObservation(*this);
}


/***********************************************************************//**
 * @brief Set LAT optimizer test methods
 ***************************************************************************/
void TestGLATOptimize::set(void)
{
    // Set test name
    name("LAT optimizers");

    // Append tests to test suite
    append(static_cast<pfunction>(&TestGLATOptimize::test_binned_optimizer_p6),
           "Test P6 binned optimizer");
    append(static_cast<pfunction>(&TestGLATOptimize::test_binned_optimizer_p7),
           "Test P7 binned optimizer");
    append(static_cast<pfunction>(&TestGLATOptimize::test_binned_optimizer_p8),
           "Test P8 binned optimizer");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGLATOptimize* TestGLATOptimize::clone(void) const
{
    // Clone test suite
    return new TestGLATOptimize(*this);
}


/***********************************************************************//**
 * @brief Test Fermi/LAT Pass 6 response handling
 ***************************************************************************/
void TestGLATResponse::test_response_p6(void)
{
    // Test Pass 6 IRFs
    test_one_response("P6_V3_DIFFUSE");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test Fermi/LAT Pass 7 response handling
 ***************************************************************************/
void TestGLATResponse::test_response_p7(void)
{
    // Test Pass 7 IRFs
    test_one_response("P7SOURCE_V6");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test Fermi/LAT Pass 8 response handling
 ***************************************************************************/
void TestGLATResponse::test_response_p8(void)
{
    // Set Pass 8 response function
    std::string irf = "P8R2_SOURCE_V6";

    // Test Pass 8 IRFs
    test_one_response(irf);

    // Test loading of Pass 8 PSF response
    test_try("Test loading of Pass 8 PSF response");
    try {
        GLATResponse rsp;
        rsp.load(irf+"::psf");
        rsp.load(irf+"::psf0");
        rsp.load(irf+"::psf1");
        rsp.load(irf+"::psf2");
        rsp.load(irf+"::psf3");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test loading of Pass 8 EDISP response
    test_try("Test loading of Pass 8 EDISP response");
    try {
        GLATResponse rsp;
        rsp.load(irf+"::edisp");
        rsp.load(irf+"::edisp0");
        rsp.load(irf+"::edisp1");
        rsp.load(irf+"::edisp2");
        rsp.load(irf+"::edisp3");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test one specific response
 *
 * @param[in] irf Instrument response function.
 *
 * Verifies the ability to load and to save Fermi/LAT response functions.
 ***************************************************************************/
void TestGLATResponse::test_one_response(const std::string& irf)
{
    // Set FITS filename
    std::string fitsfile = "test_rsp_" + irf + ".fits";
    
    // Remove FITS file
    std::string cmd = "rm -rf " + fitsfile;
    system(cmd.c_str());

    // Try loading the response
    test_try("Test loading the response");
    try {
        GLATResponse rsp;
        rsp.load(irf+"::front");
        rsp.load(irf+"::back");
        rsp.load(irf);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Try saving the response
    test_try("Test saving the response");
    try {
        GLATResponse rsp;
        rsp.load(irf);
        rsp.save(fitsfile);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test livetime cube handling
 *
 * Verifies handling of Fermi/LAT Pass 6 livetime cube.
 ***************************************************************************/
void TestGLATLtCube::test_ltcube_p6(void)
{
    // Test Pass 6 livetime cube
    test_one_ltcube(dirPass6, 339733.528629);

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test livetime cube handling
 *
 * Verifies handling of Fermi/LAT Pass 7 livetime cube.
 ***************************************************************************/
void TestGLATLtCube::test_ltcube_p7(void)
{
    // Test various datasets
    test_one_ltcube(dirPass7, 248009.734604);

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test livetime cube handling
 *
 * Verifies handling of Fermi/LAT Pass 8 livetime cube.
 ***************************************************************************/
void TestGLATLtCube::test_ltcube_p8(void)
{
    // Test various datasets
    test_one_ltcube(dirPass8, 31911.50386047);

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test livetime cube handling for a specific dataset
 *
 * @param[in] datadir Directory of test data.
 * @param[in] reference Reference value.
 *
 * Verifies the ability to handle Fermi/LAT livetime cubes.
 ***************************************************************************/
void TestGLATLtCube::test_one_ltcube(const std::string& datadir, const double& reference)
{
    // Set filenames
    std::string lat_ltcube = datadir+"/ltcube.fits";
    std::string file1      = "test_lat_ltcube.fits";
    std::string file2      = "test_lat_ltcube_phi.fits";

    // Load livetime cube
    test_try("Load livetime cube");
    try {
        // Load livetime cube
        GLATLtCube ltcube(lat_ltcube);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Load livetime cube
    GLATLtCube ltcube(lat_ltcube);

    // Initialise sky direction and energy
    GSkyDir dir;
    GEnergy energy;

    // Test cos theta integration operator
    double sum = ltcube(dir, energy, test_fct1);
    test_value(sum, reference, 0.001, "Livetime cube sum computation");

    // Test cos theta and phi integration operator. The sum differs from
    // above since the actual test dataset does not cover the same time
    // interval
    sum = ltcube(dir, energy, test_fct2);
    test_value(sum, 0.0, 0.001, "Livetime cube sum computation");

    // Create livetime skymap (no phi dependence)
    test_try("Create livetime skymap (no phi dependence)");
    try {
        GSkyMap map("GAL", 64, "RING", 1);
        GLATLtCube ltcube(lat_ltcube);
        GEnergy energy;
        for (int i = 0; i < map.npix(); ++i) {
            GSkyDir dir = map.inx2dir(i);
            map(i) = ltcube(dir, energy, test_fct1);
        }
        map.save(file1, true);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Create livetime skymap (phi dependence)
    test_try("Create livetime skymap (phi dependence)");
    try {
        GSkyMap map("GAL", 64, "RING", 1);
        GLATLtCube ltcube(lat_ltcube);
        GEnergy energy;
        for (int i = 0; i < map.npix(); ++i) {
            GSkyDir dir = map.inx2dir(i);
            map(i) = ltcube(dir, energy, test_fct2);
        }
        map.save(file2, true);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test unbinned observation handling
 *
 * Verifies handling of Pass 6 unbinned data. 
 ***************************************************************************/
void TestGLATObservation::test_unbinned_obs_p6(void)
{
    // Test various datasets
    test_one_unbinned_obs(dirPass6);

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test unbinned observation handling
 *
 * Verifies handling of Pass 7 unbinned data. 
 ***************************************************************************/
void TestGLATObservation::test_unbinned_obs_p7(void)
{
    // Test various datasets
    test_one_unbinned_obs(dirPass7);

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test unbinned observation handling
 *
 * Verifies handling of Pass 8 unbinned data. 
 ***************************************************************************/
void TestGLATObservation::test_unbinned_obs_p8(void)
{
    // Test various datasets
    test_one_unbinned_obs(dirPass8);

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test binned observation handling
 *
 * Verifies the ability to handle Pass 6 binned Fermi/LAT data.
 ***************************************************************************/
void TestGLATObservation::test_binned_obs_p6(void)
{
    // Test various datasets
    test_one_binned_obs(dirPass6, "P6_V3_DIFFUSE");

    // Exit test
    return;

}


/***********************************************************************//**
 * @brief Test binned observation handling
 *
 * Verifies the ability to handle Pass 7 binned Fermi/LAT data.
 ***************************************************************************/
void TestGLATObservation::test_binned_obs_p7(void)
{
    // Test various datasets
    test_one_binned_obs(dirPass7, "P7SOURCE_V6");

    // Exit test
    return;

}


/***********************************************************************//**
 * @brief Test binned observation handling
 *
 * Verifies the ability to handle Pass 8 binned Fermi/LAT data.
 ***************************************************************************/
void TestGLATObservation::test_binned_obs_p8(void)
{
    // Test various datasets
    test_one_binned_obs(dirPass8, "P8R2_SOURCE_V6");

    // Exit test
    return;

}


/***********************************************************************//**
 * @brief Test unbinned observation handling for a specific dataset
 *
 * @param[in] datadir Directory of test data.
 *
 * Verifies the ability to handle unbinned Fermi/LAT data.
 ***************************************************************************/
void TestGLATObservation::test_one_unbinned_obs(const std::string& datadir)
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

    // Try loading event list
    GLATEventList list(lat_ft1);
    test_value(list.number(), nevents, "Test number of events in list.");

    // Load unbinned LAT observation
    test_try("Load unbinned LAT observation");
    try {
        run.load_unbinned(lat_ft1, lat_ft2, "");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Add observation (twice) to data
    test_try("Append observation twice");
    try {
        run.id("0001");
        obs.append(run);
        run.id("0002");
        obs.append(run);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Loop over all events
    const GEvents *ptr = run.events();
    int num = 0;
    for (int i = 0; i < ptr->size(); ++i) {
        num++;
    }
    test_value(num, nevents, 1.0e-20, "Test event iterator");

    // Test XML loading
    test_try("Test XML loading");
    try {
        obs = GObservations(lat_unbin_xml);
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
 * @brief Test binned observation handling for a specific dataset
 *
 * @param[in] datadir Directory of test data.
 * @param[in] irf Instrument response function.
 *
 * Verifies the ability to handle binned Fermi/LAT data.
 ***************************************************************************/
void TestGLATObservation::test_one_binned_obs(const std::string& datadir,
                                              const std::string& irf)
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

    // Try loading event list
    GLATEventCube cube(lat_cntmap);
    test_value(double(cube.number()), nevents, "Check number of events in cube.");

    // Load LAT binned observation from counts map
    test_try("Load LAT binned observation");
    try {
        run.load_binned(lat_cntmap, "", "");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Reload LAT binned observation from source map
    test_try("Reload LAT binned observation");
    try {
        run.load_binned(lat_srcmap, "", "");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Add observation (twice) to data
    test_try("Append observation twice");
    try {
        run.id("0001");
        obs.append(run);
        run.id("0002");
        obs.append(run);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Loop over all events using iterator
    const GEvents* events = run.events();
    int num = 0;
    int sum = 0;
    for (int i = 0; i < events->size(); ++i) {
        num++;
        sum += (int)((*events)[i]->counts());
    }
    test_value(sum, nevents, 1.0e-20, "Test event iterator (counts)");
    test_value(num, nsize, 1.0e-20, "Test event iterator (bins)");

    // Test mean PSF
    test_try("Test mean PSF");
    try {
        run.load_binned(lat_srcmap, lat_expmap, lat_ltcube);
        run.response(irf);
        GSkyDir dir;
        GLATMeanPsf psf(dir, run);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test XML loading
    test_try("Test XML loading");
    try {
        obs = GObservations(lat_bin_xml);
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
 * @brief Test binned optimizer handling
 *
 * Verifies the ability to handle binned Pass 6 Fermi/LAT optimization.
 ***************************************************************************/
void TestGLATOptimize::test_binned_optimizer_p6(void)
{
    // Set expected fit results
    double fit_results[] = {1, 0,
                            2.34236, 0.7339917,
                            1, 0,
                            1, 0,
                            0.8167256818, 0.08569407542,
                            1, 0,
                            83.6331, 0,
                            22.0145, 0,
                            2.103749294e-06, 1.560611406e-07,
                            -2.194252857, 0.05165097,
                            100, 0,
                            500000, 0,
                            1, 0};

    // Test various datasets
    test_one_binned_optimizer(dirPass6, "P6_V3_DIFFUSE", fit_results);

    // Exit test
    return;

}


/***********************************************************************//**
 * @brief Test binned optimizer handling
 *
 * Verifies the ability to handle binned Pass 7 Fermi/LAT optimization.
 ***************************************************************************/
void TestGLATOptimize::test_binned_optimizer_p7(void)
{
    // Set expected fit results
    double fit_results[] = {1, 0,
                            2.37468, 0.4548979,
                            1, 0,
                            1, 0,
                            0.8419824722, 0.06340531926,
                            1, 0,
                            83.6331, 0,
                            22.0145, 0,
                            1.922525774e-06, 1.209507237e-07,
                            -2.12421, 0.0493105484,
                            100, 0,
                            500000, 0,
                            1, 0};

    // Test various datasets
    test_one_binned_optimizer(dirPass7, "P7SOURCE_V6", fit_results);

    // Exit test
    return;

}


/***********************************************************************//**
 * @brief Test binned optimizer handling
 *
 * Verifies the ability to handle binned Pass 8 Fermi/LAT optimization.
 ***************************************************************************/
void TestGLATOptimize::test_binned_optimizer_p8(void)
{
    // Set expected fit results
    double fit_results[] = {1, 0,
//                            2.37468, 0.4548979, // Pass 7
                            1.99669, 1.00881,   // Pass 8
                            1, 0,
                            1, 0,
//                            0.8419824722, 0.06340531926, // Pass 7
                            1.00878, 0.0940645,          // Pass 8
                            1, 0,
                            83.6331, 0,
                            22.0145, 0,
                            1.922525774e-06, 1.209507237e-07,
//                            -2.12421, 0.0493105484, // Pass 7
                            -2.31512, 0.147089,     // Pass 8
                            100, 0,
                            500000, 0,
                            1, 0};

    // Test various datasets
    test_one_binned_optimizer(dirPass8, "P8R2_SOURCE_V6", fit_results);

    // Exit test
    return;

}


/***********************************************************************//**
 * @brief Test binned optimizer
 *
 * @param[in] datadir Directory of test data.
 * @param[in] irf Instrument response function.
 * @param[in] fit_results Expected fit result.
 *
 * Verifies the ability optimize binned Fermi/LAT data.
 ***************************************************************************/
void TestGLATOptimize::test_one_binned_optimizer(const std::string& datadir,
                                                 const std::string& irf,
                                                 const double*      fit_results)
{
    // Set filenames
    std::string lat_srcmap    = datadir+"/srcmap.fits";
    std::string lat_expmap    = datadir+"/binned_expmap.fits";
    std::string lat_ltcube    = datadir+"/ltcube.fits";
    std::string lat_model_xml = datadir+"/source_model.xml";

    // Setup GObservations for optimizing
    GObservations   obs;
    GLATObservation run;
    test_try("Setup for optimization");
    try {
        run.load_binned(lat_srcmap, lat_expmap, lat_ltcube);
        run.response(irf);
        obs.append(run);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Load models from XML file
    obs.models(lat_model_xml);

    // Setup LM optimizer
    test_try("Perform LM optimization");
    try {
        GOptimizerLM opt;
        opt.max_iter(1000);
        obs.optimize(opt);
        obs.errors(opt);
        test_try_success();
        for (int i = 0, j = 0; i < obs.models().size(); ++i) {
            const GModel* model = obs.models()[i];
            for (int k = 0; k < model->size(); ++k) {
                GModelPar par  = (*model)[k];
                std::string msg = "Verify optimization result for " + par.print();
                test_value(par.value(), fit_results[j++], 5.0e-5, msg);
                test_value(par.error(), fit_results[j++], 5.0e-5, msg);
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
    GTestSuites testsuites("LAT instrument specific class testing");

    // Check if data directory exists
    bool has_data = (access(datadir.c_str(), R_OK) == 0);

    // Set CALDB environment variable
    setenv("CALDB", lat_caldb.c_str(), 1);

    // Initially assume that we pass all tests
    bool success = true;

    // Create test suites and append them to the container
    TestGLATResponse    rsp;
    TestGLATLtCube      ltcube;
    TestGLATObservation obs;
    TestGLATOptimize    opt;
    testsuites.append(rsp);
    if (has_data) {
        testsuites.append(ltcube);
        testsuites.append(obs);
        testsuites.append(opt);
    }

    // Run the testsuites
    success = testsuites.run();

    // Save test report
    testsuites.save("reports/GLAT.xml");

    // Return success status
    return (success ? 0 : 1);
}
