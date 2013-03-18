/***************************************************************************
 *                       test_COM.cpp - test COM classes                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2013 by Juergen Knoedlseder                         *
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
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <unistd.h>
#include "GTools.hpp"
#include "test_COM.hpp"

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
 * @brief Set COMPTEL response test methods
 ***************************************************************************/
void TestGCOMResponse::set(void)
{
    // Set test name
    name("GCOMResponse");

    // Append tests to test suite
    append(static_cast<pfunction>(&TestGCOMResponse::test_inst_dir), "Test instrument direction");
    append(static_cast<pfunction>(&TestGCOMResponse::test_pointing), "Test pointing");
    append(static_cast<pfunction>(&TestGCOMResponse::test_response), "Test response");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set COMPTEL observation test methods
 ***************************************************************************/
void TestGCOMObservation::set(void)
{
    // Set test name
    name("GCOMObservation");

    // Append tests to test suite
    append(static_cast<pfunction>(&TestGCOMObservation::test_binned_obs), "Test binned observation");
    append(static_cast<pfunction>(&TestGCOMObservation::test_event_bin), "Test event bin");
    append(static_cast<pfunction>(&TestGCOMObservation::test_event_cube), "Test event cube");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set COMPTEL optimizer test methods
 ***************************************************************************/
void TestGCOMOptimize::set(void)
{
    // Set test name
    name("COMPTEL optimizers");

    // Append tests to test suite
    append(static_cast<pfunction>(&TestGCOMOptimize::test_binned_optimizer), "Test binned optimizer");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GCOMInstDir class
 ***************************************************************************/
void TestGCOMResponse::test_inst_dir(void)
{
    // Test constructors
    test_try("Test constructors");
    try {
        // Void constructor
        GCOMInstDir dir1;

        // Copy constructor
        GCOMInstDir dir2(dir1);

        // If we arrived here, signal success
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Create void object
    GCOMInstDir dir;

    // skydir method
    GSkyDir sky;
    sky.radec_deg(37.0, 45.3);
    dir.dir(sky);
    test_assert(dir.dir() == sky, "Test skydir() method.",
                "Expected "+sky.print()+", found "+dir.dir().print());

    // phibar method
    dir.phibar(27.2);
    test_value(dir.phibar(), 27.2, 1.0e-10, "Test phibar() method.");

    // copy constructor
    GCOMInstDir dir_copy(dir);
    test_assert(dir_copy.dir() == sky, "Test copy constructor method.",
                "Expected "+sky.print()+", found "+dir_copy.dir().print());
    test_value(dir_copy.phibar(), 27.2, 1.0e-10, "Test copy constructor method.");

    // assignment operator
    GCOMInstDir dir_assign = dir;
    test_assert(dir_assign.dir() == sky, "Test assignment operator method.",
                "Expected "+sky.print()+", found "+dir_assign.dir().print());
    test_value(dir_assign.phibar(), 27.2, 1.0e-10, "Test assignment operator method.");

    // clone
    GCOMInstDir* dir_clone = dir.clone();
    test_assert(dir_clone->dir() == sky, "Test clone() method.",
                "Expected "+sky.print()+", found "+dir_clone->dir().print());
    test_value(dir_clone->phibar(), 27.2, 1.0e-10, "Test clone() method.");

    // clear
    dir.clear();
    sky.clear();
    test_assert(dir.dir() == sky, "Test clear() method.",
                "Expected "+sky.print()+", found "+dir.dir().print());
    test_value(dir.phibar(), 0.0, 1.0e-10, "Test clean() method.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GCOMPointing class
 ***************************************************************************/
void TestGCOMResponse::test_pointing(void)
{
    // Test constructors
    test_try("Test constructors");
    try {
        // Void constructor
        GCOMPointing pnt1;

        // Copy constructor
        GCOMPointing pnt2(pnt1);

        // If we arrived here, signal success
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Create void object
    GCOMPointing pnt;

    // dir method
    GSkyDir sky;
    sky.radec_deg(37.0, 45.3);
    pnt.dir(sky);
    test_assert(pnt.dir() == sky, "Test dir() method.",
                "Expected "+sky.print()+", found "+pnt.dir().print());

    // copy constructor
    GCOMPointing pnt_copy(pnt);
    test_assert(pnt_copy.dir() == sky, "Test copy constructor method.",
                "Expected "+sky.print()+", found "+pnt_copy.dir().print());

    // assignment operator
    GCOMPointing pnt_assign = pnt;
    test_assert(pnt_assign.dir() == sky, "Test assignment operator method.",
                "Expected "+sky.print()+", found "+pnt_assign.dir().print());

    // clone
    GCOMPointing* pnt_clone = pnt.clone();
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
void TestGCOMResponse::test_response(void)
{
    // Test constructors
    test_try("Test constructors");
    try {
        // Void constructor
        GCOMResponse rsp1;

        // Copy constructor
        GCOMResponse rsp2(rsp1);

        // If we arrived here, signal success
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test response loading
    test_try("Test response loading");
    try {
        // Construct response from datasets
        GCOMResponse rsp(com_iaq, com_caldb);

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
void TestGCOMObservation::test_binned_obs(void)
{
    // Test constructors
    test_try("Test constructors");
    try {
        // Void constructor
        GCOMObservation obs1;

        // Copy constructor
        GCOMObservation obs2(obs1);

        // If we arrived here, signal success
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test observation loading
    test_try("Test observation loading");
    try {
        // Construct observation from datasets
        GCOMObservation obs(com_dre, com_drb, com_drg, com_drx);
//std::cout << obs << std::endl;

        // If we arrived here, signal success
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test XML constructor
    test_try("Test XML constructor");
    try {
        // Construct observation from datasets
        GObservations obs(com_obs);
//std::cout << obs << std::endl;

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
 * @brief Checks handling of COMPTEL event bin
 ***************************************************************************/
void TestGCOMObservation::test_event_bin(void)
{
    // Test COMPTEL event bin methods (one by one)
    test_try("Test event bin methods");
    try {
        // Event bin void constructor
        GCOMEventBin bin;

        // Copy constructor
        GCOMEventBin bin2(bin);

        // Assignment operator
        GCOMEventBin bin3 = bin;

        // clear method
        bin.clear();

        // clone method
        GCOMEventBin* bin4 = bin.clone();

        // size method
        test_value(bin.size(), 0.0, 1.0e-10, "Test size() method.");
        test_value(bin4->size(), 0.0, 1.0e-10, "Test size() method.");

        // dir method
        GCOMInstDir dir = bin.dir();

        // energy method
        GEnergy energy = bin.energy();

        // time method
        GTime time = bin.time();

        // counts method
        test_value(bin.counts(), 0.0, 1.0e-10, "Test counts() method for zero counts.");

        // error method
        test_value(bin.error(), 0.0, 1.0e-10, "Test error() method for zero counts.");

        // counts setting
        bin.counts(1.0);
        test_value(bin.counts(), 1.0, 1.0e-10, "Test counts() method for one count.");
        test_value(bin.error(), 1.0, 1.0e-10, "Test error() method for one count.");

        // print method
        std::string text = bin.print();
        test_assert(text == "1", "Test print() method.");

        // omega method
        test_value(bin.omega(), 0.0, 1.0e-10, "Test omega() method.");

        // ewidth method
        GEnergy ewidth = bin.ewidth();

        // ontime method
        test_value(bin.ontime(), 0.0, 1.0e-10, "Test ontime() method.");

        // If we arrived here, signal success
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test COMPTEL event bin operations
    GCOMEventBin bin;
    bin.counts(3.3);
    test_value(bin.counts(), 3.3, 1.0e-10, "Test counts() method for 3.3 counts.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Checks handling of COMPTEL event cube
 ***************************************************************************/
void TestGCOMObservation::test_event_cube(void)
{
    // Event cube void constructor
    GCOMEventCube cube;

    // Event cube load constructor
    GCOMEventCube cube2(com_dre);

    // Event cube copy constructor
    GCOMEventCube cube3(cube2);
    test_value(cube2.size(), cube3.size(), "Test copy constructor.",
               "Different cube size after using the"
               " copy constructor (before="+str(cube2.size())+
               " after="+str(cube3.size())+")");
    test_value(cube2.number(), cube3.number(), "Test copy constructor.",
               "Different number of events after using the"
               " copy constructor (before="+str(cube2.number())+
               " after="+str(cube3.number())+").");

    // Event cube assignment operator
    GCOMEventCube cube4 = cube2;
    test_value(cube2.size(), cube4.size(), "Test assignment operator.",
               "Different cube size after using the"
               " assignment operator (before="+str(cube2.size())+
               " after="+str(cube4.size())+")");
    test_value(cube2.number(), cube4.number(), "Test assignment operator.",
               "Different number of events after using the"
               " assignment operator (before="+str(cube2.number())+
               " after="+str(cube4.number())+").");
        
    // clear method
    cube4.clear();
    test_value(cube4.size(), 0, "Test clear method.",
               "Expected event cube with size 0 after"
               " clear but found size "+str(cube4.size())+".");
    test_value(cube4.number(), 0, "Test clear method.",
               "Expected 0 events in cube after"
               " clear but found "+str(cube4.size())+" events.");

    // clone method
    GCOMEventCube* cube5 = cube2.clone();
    test_value(cube2.size(), cube5->size(), "Test clone() method.",
               "Different cube size after cloning"
               " (before="+str(cube2.size())+
               " after="+str(cube5->size())+")");
    test_value(cube2.number(), cube5->number(), "Test clone() method.",
               "Different number of events after cloning"
               " (before="+str(cube2.number())+
               " after="+str(cube5->number())+").");

    // size method
    test_value(cube2.size(), 140600, "Test size() method.",
               "Expected cube dimension 140600, found "+
               str(cube2.size())+".");

    // dim method
    test_value(cube2.dim(), 3, "Test dim() method.",
               "Expected 3 cube dimensions, found "+
               str(cube2.dim())+".");

    // naxis method
    test_value(cube2.naxis(0), 76, "Test naxis(0) method.",
               "Expected Chi axis dimension 76, found "+
               str(cube2.naxis(0))+".");
    test_value(cube2.naxis(1), 74, "Test naxis(1) method.",
               "Expected Chi axis dimension 74, found "+
               str(cube2.naxis(1))+".");
    test_value(cube2.naxis(2), 25, "Test naxis(2) method.",
               "Expected Chi axis dimension 25, found "+
               str(cube2.naxis(2))+".");

    // nchi, npsi, nphi methods
    test_value(cube2.nchi(), 76, "Test nchi() method.",
               "Expected Chi axis dimension 76, found "+
               str(cube2.nchi())+".");
    test_value(cube2.npsi(), 74, "Test npsi() method.",
               "Expected Chi axis dimension 74, found "+
               str(cube2.npsi())+".");
    test_value(cube2.nphi(), 25, "Test nphi() method.",
               "Expected Chi axis dimension 25, found "+
               str(cube2.nphi())+".");

    // npix method
    int npix = cube2.nchi() * cube2.npsi();
    test_value(cube2.npix(), npix, "Test npix() method.",
               "Expected "+str(npix)+" pixels in (Chi,Psi) plane, found "+
               str(cube2.npix())+".");

    // number method
    test_value(cube2.number(), 316141, "Test number() method.",
               "Expected 316141 events in cube, found "+
               str(cube2.number())+".");

    // event access
    double sum = 0.0;
    for (int i = 0; i < cube2.size(); ++i) {
        sum += cube2[i]->counts();
    }
    test_value(cube2.number(), int(sum+0.5), "Test event access.",
               "Expected "+str(cube2.number())+
               " events in cube, found "+
               str(int(sum+0.5))+" by summing over all elements.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test binned optimizer
 ***************************************************************************/
void TestGCOMOptimize::test_binned_optimizer(void)
{
    // Declare observations
    GObservations   obs;
    GCOMObservation run;

    // Load binned COMPTEL observation
    test_try("Load binned COMPTEL observation");
    try {
        run.load(com_dre, com_drb, com_drg, com_drx);
        run.response(com_iaq, com_caldb);
        obs.append(run);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Load models from XML file
    obs.models(com_model);

    // Perform LM optimization
    double fit_results[] = {83.4258, 0.157307,
                            21.5953, 0.143873,
                            0.001907, 8.31055e-05,
                            -2.05, 0,
                            1, 0,
                            1, 0,
                            1, 0,
                            0, 0,
                            3, 0,
                            0, 0,
                            5, 0,
                            3.15461, 0.0904736,
                            7, 0,
                            27.8763, 0.272976,
                            9, 0,
                            43.6555, 0.348297,
                            11, 0,
                            54.1233, 0.396509,
                            13, 0,
                            60.9322, 0.431667,
                            15, 0,
                            65.1235, 0.457003,
                            17, 0,
                            64.8946, 0.466807,
                            19, 0,
                            60.7385, 0.458565,
                            21, 0,
                            57.5857, 0.451466,
                            23, 0,
                            54.447,  0.441844,
                            25, 0,
                            51.9042, 0.433244,
                            27, 0,
                            49.0977, 0.422531,
                            29, 0,
                            47.1811, 0.416003,
                            31, 0,
                            44.2415, 0.406013,
                            33, 0,
                            43.1503, 0.40499,
                            35, 0,
                            41.7076, 0.403297,
                            37, 0,
                            39.92,   0.39994,
                            39, 0,
                            37.602,  0.39392,
                            41, 0,
                            36.7116, 0.395238,
                            43, 0,
                            35.8177, 0.397315,
                            45, 0,
                            33.7205, 0.392628,
                            47, 0,
                            29.9617, 0.377506,
                            49, 0,
                            27.9992, 0.372667};
    test_try("Perform LM optimization");
    try {
        GOptimizerLM opt;
        opt.max_iter(100);
        obs.optimize(opt);
        test_try_success();
        for (int i = 0, j = 0; i < obs.models().size(); ++i) {
            const GModel* model = obs.models()[i];
            for (int k = 0; k < model->size(); ++k) {
                GModelPar par  = (*model)[k];
                std::string msg = "Verify optimization result for " + par.print();
                test_value(par.Value(), fit_results[j++], 5.0e-5, msg);
                test_value(par.Error(), fit_results[j++], 5.0e-5, msg);
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
    GTestSuites testsuites("COMPTEL instrument specific class testing");

    // Set GAMMALIB_CALDB environment variable
    std::string caldb = "GAMMALIB_CALDB="+com_caldb;
    putenv((char*)caldb.c_str());

    // Initially assume that we pass all tests
    bool success = true;

    // Create test suites and append them to the container
    TestGCOMResponse    rsp;
    TestGCOMObservation obs;
    TestGCOMOptimize    opt;
    testsuites.append(rsp);
    testsuites.append(obs);
    testsuites.append(opt);

    // Run the testsuites
    success = testsuites.run();

    // Save test report
    testsuites.save("reports/GCOM.xml");

    // Return success status
    return (success ? 0 : 1);
}
