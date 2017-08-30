/***************************************************************************
 *                       test_COM.cpp - test COM classes                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2017 by Juergen Knoedlseder                         *
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
#include <cstdlib>     // getenv
#include "GTools.hpp"
#include "test_COM.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Constants __________________________________________________________ */
const std::string datadir   = std::getenv("TEST_COM_DATA");
const std::string com_caldb = datadir + "/../../caldb";
const std::string com_iaq   = "ENERG(1.0-3.0)MeV";        // 1-3 MeV
const std::string com_dre   = datadir+"/m50439_dre.fits"; // 1-3 MeV
const std::string com_drb   = datadir+"/m34997_drg.fits";
const std::string com_drg   = datadir+"/m34997_drg.fits";
const std::string com_drx   = datadir+"/m32171_drx.fits";
const std::string com_obs   = datadir+"/obs.xml";
const std::string com_model = datadir+"/crab.xml";


/***********************************************************************//**
 * @brief Set COMPTEL response test methods
 ***************************************************************************/
void TestGCOMResponse::set(void)
{
    // Set test name
    name("GCOMResponse");

    // Append tests to test suite
    append(static_cast<pfunction>(&TestGCOMResponse::test_inst_dir),
           "Test instrument direction");
    append(static_cast<pfunction>(&TestGCOMResponse::test_response),
           "Test response");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGCOMResponse* TestGCOMResponse::clone(void) const
{
    // Clone test suite
    return new TestGCOMResponse(*this);
}


/***********************************************************************//**
 * @brief Set COMPTEL observation test methods
 ***************************************************************************/
void TestGCOMObservation::set(void)
{
    // Set test name
    name("GCOMObservation");

    // Append tests to test suite
    append(static_cast<pfunction>(&TestGCOMObservation::test_binned_obs),
           "Test binned observation");
    append(static_cast<pfunction>(&TestGCOMObservation::test_event_bin),
           "Test event bin");
    append(static_cast<pfunction>(&TestGCOMObservation::test_event_cube),
           "Test event cube");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGCOMObservation* TestGCOMObservation::clone(void) const
{
    // Clone test suite
    return new TestGCOMObservation(*this);
}


/***********************************************************************//**
 * @brief Set COMPTEL optimizer test methods
 ***************************************************************************/
void TestGCOMOptimize::set(void)
{
    // Set test name
    name("COMPTEL optimizers");

    // Append tests to test suite
    append(static_cast<pfunction>(&TestGCOMOptimize::test_binned_optimizer),
           "Test binned optimizer");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGCOMOptimize* TestGCOMOptimize::clone(void) const
{
    // Clone test suite
    return new TestGCOMOptimize(*this);
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
    delete dir_clone;

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
        GCOMResponse rsp(GCaldb("cgro", "comptel"), com_iaq);

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
    // Test event bin void constructor
    GCOMEventBin bin;
    test_value(bin.classname(), "GCOMEventBin", "Check classname() for empty bin");
    test_value(bin.size(), 0.0, 1.0e-10, "Check size() for empty bin");
    test_value(bin.dir().dir().ra_deg(), 0.0, 1.0e-10,
               "Check dir().dir().ra_deg() for empty bin");
    test_value(bin.dir().dir().dec_deg(), 0.0, 1.0e-10,
               "Check dir().dir().dec_deg() for empty bin");
    test_value(bin.dir().phibar(), 0.0, 1.0e-10,
               "Check dir().phibar() for empty bin");
    test_value(bin.energy().MeV(), 0.0, 1.0e-10, "Check energy() for empty bin");
    test_value(bin.time().secs(), 0.0, 1.0e-10, "Check time() for empty bin");
    test_value(bin.counts(), 0.0, 1.0e-10, "Check counts() for empty bin");
    test_value(bin.error(),  0.0, 1.0e-10, "Check error() for empty bin");
    test_value(bin.index(), -1, "Check index() for empty bin");
    test_value(bin.solidangle(), 0.0, 1.0e-10, "Check solidangle() for empty bin");
    test_value(bin.ewidth().MeV(), 0.0, 1.0e-10, "Check ewidth() for empty bin");
    test_value(bin.ontime(), 0.0, 1.0e-10, "Check ontime() for empty bin");
    test_value(bin.print(), "0", "Check print() for empty bin");

    // Test setting of bin attributes (for the moment we can only set the
    // number of counts, but this should change in the future)
    bin.counts(4.0);
    test_value(bin.size(), 0.0, 1.0e-10, "Check size() for filled bin");
    test_value(bin.dir().dir().ra_deg(), 0.0, 1.0e-10,
               "Check dir().dir().ra_deg() for filled bin");
    test_value(bin.dir().dir().dec_deg(), 0.0, 1.0e-10,
               "Check dir().dir().dec_deg() for filled bin");
    test_value(bin.dir().phibar(), 0.0, 1.0e-10,
               "Check dir().phibar() for filled bin");
    test_value(bin.energy().MeV(), 0.0, 1.0e-10,
               "Check energy() for filled bin");
    test_value(bin.time().secs(), 0.0, 1.0e-10, "Check time() for filled bin");
    test_value(bin.counts(), 4.0, 1.0e-10, "Check counts() for filled bin");
    test_value(bin.error(),  2.0, 1.0e-10, "Check error() for filled bin");
    test_value(bin.index(), -1, "Check index() for filled bin");
    test_value(bin.solidangle(), 0.0, 1.0e-10,
               "Check solidangle() for filled bin");
    test_value(bin.ewidth().MeV(), 0.0, 1.0e-10,
               "Check ewidth() for filled bin");
    test_value(bin.ontime(), 0.0, 1.0e-10, "Check ontime() for filled bin");
    test_value(bin.print(), "4", "Check print() for filled bin");

    // Test copy constructor
    GCOMEventBin bin2(bin);
    test_value(bin2.size(), 0.0, 1.0e-10, "Check size() for copied bin");
    test_value(bin2.dir().dir().ra_deg(), 0.0, 1.0e-10,
               "Check dir().dir().ra_deg() for copied bin");
    test_value(bin2.dir().dir().dec_deg(), 0.0, 1.0e-10,
               "Check dir().dir().dec_deg() for copied bin");
    test_value(bin2.dir().phibar(), 0.0, 1.0e-10,
               "Check dir().phibar() for copied bin");
    test_value(bin2.energy().MeV(), 0.0, 1.0e-10,
               "Check energy() for copied bin");
    test_value(bin2.time().secs(), 0.0, 1.0e-10, "Check time() for copied bin");
    test_value(bin2.counts(), 4.0, 1.0e-10, "Check counts() for copied bin");
    test_value(bin2.error(),  2.0, 1.0e-10, "Check error() for copied bin");
    test_value(bin2.index(), -1, "Check index() for copied bin");
    test_value(bin2.solidangle(), 0.0, 1.0e-10,
               "Check solidangle() for copied bin");
    test_value(bin2.ewidth().MeV(), 0.0, 1.0e-10,
               "Check ewidth() for copied bin");
    test_value(bin2.ontime(), 0.0, 1.0e-10, "Check ontime() for copied bin");
    test_value(bin2.print(), "4", "Check print() for copied bin");

    // Assignment operator
    GCOMEventBin bin3 = bin;
    test_value(bin3.size(), 0.0, 1.0e-10, "Check size() for assigned bin");
    test_value(bin3.dir().dir().ra_deg(), 0.0, 1.0e-10,
               "Check dir().dir().ra_deg() for assigned bin");
    test_value(bin3.dir().dir().dec_deg(), 0.0, 1.0e-10,
               "Check dir().dir().dec_deg() for assigned bin");
    test_value(bin3.dir().phibar(), 0.0, 1.0e-10,
               "Check dir().phibar() for assigned bin");
    test_value(bin3.energy().MeV(), 0.0, 1.0e-10,
               "Check energy() for assigned bin");
    test_value(bin3.time().secs(), 0.0, 1.0e-10,
               "Check time() for assigned bin");
    test_value(bin3.counts(), 4.0, 1.0e-10, "Check counts() for assigned bin");
    test_value(bin3.error(),  2.0, 1.0e-10, "Check error() for assigned bin");
    test_value(bin3.index(), -1, "Check index() for assigned bin");
    test_value(bin3.solidangle(), 0.0, 1.0e-10,
               "Check solidangle() for assigned bin");
    test_value(bin3.ewidth().MeV(), 0.0, 1.0e-10,
               "Check ewidth() for assigned bin");
    test_value(bin3.ontime(), 0.0, 1.0e-10, "Check ontime() for assigned bin");
    test_value(bin3.print(), "4", "Check print() for assigned bin");

    // clone method
    GCOMEventBin* bin4 = bin.clone();
    test_value(bin4->size(), 0.0, 1.0e-10, "Check size() for cloned bin");
    test_value(bin4->dir().dir().ra_deg(), 0.0, 1.0e-10,
               "Check dir().dir().ra_deg() for cloned bin");
    test_value(bin4->dir().dir().dec_deg(), 0.0, 1.0e-10,
               "Check dir().dir().dec_deg() for cloned bin");
    test_value(bin4->dir().phibar(), 0.0, 1.0e-10,
               "Check dir().phibar() for cloned bin");
    test_value(bin4->energy().MeV(), 0.0, 1.0e-10,
               "Check energy() for cloned bin");
    test_value(bin4->time().secs(), 0.0, 1.0e-10,
               "Check time() for cloned bin");
    test_value(bin4->counts(), 4.0, 1.0e-10, "Check counts() for cloned bin");
    test_value(bin4->error(),  2.0, 1.0e-10, "Check error() for cloned bin");
    test_value(bin4->index(), -1, "Check index() for cloned bin");
    test_value(bin4->solidangle(), 0.0, 1.0e-10,
               "Check solidangle() for cloned bin");
    test_value(bin4->ewidth().MeV(), 0.0, 1.0e-10,
               "Check ewidth() for cloned bin");
    test_value(bin4->ontime(), 0.0, 1.0e-10, "Check ontime() for cloned bin");
    test_value(bin4->print(), "4", "Check print() for cloned bin");

    // clear method
    bin.clear();
    test_value(bin.classname(), "GCOMEventBin",
               "Check classname() for cleared bin");
    test_value(bin.size(), 0.0, 1.0e-10, "Check size() for cleared bin");
    test_value(bin.dir().dir().ra_deg(), 0.0, 1.0e-10,
               "Check dir().dir().ra_deg() for cleared bin");
    test_value(bin.dir().dir().dec_deg(), 0.0, 1.0e-10,
               "Check dir().dir().dec_deg() for cleared bin");
    test_value(bin.dir().phibar(), 0.0, 1.0e-10,
               "Check dir().phibar() for cleared bin");
    test_value(bin.energy().MeV(), 0.0, 1.0e-10,
               "Check energy() for cleared bin");
    test_value(bin.time().secs(), 0.0, 1.0e-10, "Check time() for cleared bin");
    test_value(bin.counts(), 0.0, 1.0e-10, "Check counts() for cleared bin");
    test_value(bin.error(),  0.0, 1.0e-10, "Check error() for cleared bin");
    test_value(bin.index(), -1, "Check index() for cleared bin");
    test_value(bin.solidangle(), 0.0, 1.0e-10,
               "Check solidangle() for cleared bin");
    test_value(bin.ewidth().MeV(), 0.0, 1.0e-10,
               "Check ewidth() for cleared bin");
    test_value(bin.ontime(), 0.0, 1.0e-10, "Check ontime() for cleared bin");
    test_value(bin.print(), "0", "Check print() for cleared bin");

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
               " copy constructor (before="+gammalib::str(cube2.size())+
               " after="+gammalib::str(cube3.size())+")");
    test_value(cube2.number(), cube3.number(), "Test copy constructor.",
               "Different number of events after using the"
               " copy constructor (before="+gammalib::str(cube2.number())+
               " after="+gammalib::str(cube3.number())+").");

    // Event cube assignment operator
    GCOMEventCube cube4 = cube2;
    test_value(cube2.size(), cube4.size(), "Test assignment operator.",
               "Different cube size after using the"
               " assignment operator (before="+gammalib::str(cube2.size())+
               " after="+gammalib::str(cube4.size())+")");
    test_value(cube2.number(), cube4.number(), "Test assignment operator.",
               "Different number of events after using the"
               " assignment operator (before="+gammalib::str(cube2.number())+
               " after="+gammalib::str(cube4.number())+").");
        
    // clear method
    cube4.clear();
    test_value(cube4.size(), 0, "Test clear method.",
               "Expected event cube with size 0 after"
               " clear but found size "+gammalib::str(cube4.size())+".");
    test_value(cube4.number(), 0, "Test clear method.",
               "Expected 0 events in cube after"
               " clear but found "+gammalib::str(cube4.size())+" events.");

    // clone method
    GCOMEventCube* cube5 = cube2.clone();
    test_value(cube2.size(), cube5->size(), "Test clone() method.",
               "Different cube size after cloning"
               " (before="+gammalib::str(cube2.size())+
               " after="+gammalib::str(cube5->size())+")");
    test_value(cube2.number(), cube5->number(), "Test clone() method.",
               "Different number of events after cloning"
               " (before="+gammalib::str(cube2.number())+
               " after="+gammalib::str(cube5->number())+").");
    delete cube5;

    // size method
    test_value(cube2.size(), 140600, "Test size() method.",
               "Expected cube dimension 140600, found "+
               gammalib::str(cube2.size())+".");

    // dim method
    test_value(cube2.dim(), 3, "Test dim() method.",
               "Expected 3 cube dimensions, found "+
               gammalib::str(cube2.dim())+".");

    // naxis method
    test_value(cube2.naxis(0), 76, "Test naxis(0) method.",
               "Expected Chi axis dimension 76, found "+
               gammalib::str(cube2.naxis(0))+".");
    test_value(cube2.naxis(1), 74, "Test naxis(1) method.",
               "Expected Chi axis dimension 74, found "+
               gammalib::str(cube2.naxis(1))+".");
    test_value(cube2.naxis(2), 25, "Test naxis(2) method.",
               "Expected Chi axis dimension 25, found "+
               gammalib::str(cube2.naxis(2))+".");

    // nchi, npsi, nphi methods
    test_value(cube2.nchi(), 76, "Test nchi() method.",
               "Expected Chi axis dimension 76, found "+
               gammalib::str(cube2.nchi())+".");
    test_value(cube2.npsi(), 74, "Test npsi() method.",
               "Expected Chi axis dimension 74, found "+
               gammalib::str(cube2.npsi())+".");
    test_value(cube2.nphi(), 25, "Test nphi() method.",
               "Expected Chi axis dimension 25, found "+
               gammalib::str(cube2.nphi())+".");

    // npix method
    int npix = cube2.nchi() * cube2.npsi();
    test_value(cube2.npix(), npix, "Test npix() method.",
               "Expected "+gammalib::str(npix)+" pixels in (Chi,Psi) plane, found "+
               gammalib::str(cube2.npix())+".");

    // number method
    test_value(cube2.number(), 316141, "Test number() method.",
               "Expected 316141 events in cube, found "+
               gammalib::str(cube2.number())+".");

    // event access
    double sum = 0.0;
    for (int i = 0; i < cube2.size(); ++i) {
        sum += cube2[i]->counts();
    }
    test_value(cube2.number(), int(sum+0.5), "Test event access.",
               "Expected "+gammalib::str(cube2.number())+
               " events in cube, found "+
               gammalib::str(int(sum+0.5))+" by summing over all elements.");

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
        run.response(GCaldb("cgro", "comptel"), com_iaq);
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
                            21.5952, 0.143873,
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
                            27.9992, 0.372603};
    test_try("Perform LM optimization");
    try {
        GOptimizerLM opt;
        opt.max_iter(100);
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
    GTestSuites testsuites("COMPTEL instrument specific class testing");

    // Set CALDB environment variable
    std::string caldb = "CALDB="+com_caldb;
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
