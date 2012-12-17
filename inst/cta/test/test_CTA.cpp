/***************************************************************************
 *                       test_CTA.cpp - Test CTA classes                   *
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
 * @file test_CTA.cpp
 * @brief Implementation of CTA classes unit tests
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <iostream>
#include <unistd.h>
#include "GCTALib.hpp"
#include "GTools.hpp"
#include "test_CTA.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Constants __________________________________________________________ */
const std::string datadir       = "../inst/cta/test/data";
const std::string cta_caldb     = "../inst/cta/caldb";
const std::string cta_irf       = "kb_E_50h_v3";
const std::string cta_events    = datadir+"/crab_events.fits.gz";
const std::string cta_cntmap    = datadir+"/crab_cntmap.fits.gz";
const std::string cta_bin_xml   = datadir+"/obs_binned.xml";
const std::string cta_unbin_xml = datadir+"/obs_unbinned.xml";
const std::string cta_model_xml = datadir+"/crab.xml";
const std::string cta_rsp_xml   = datadir+"/rsp_models.xml";


/***********************************************************************//**
 * @brief Set CTA response test methods
 ***************************************************************************/
void TestGCTAResponse::set(void)
{
    // Set test name
    name("GCTAResponse");

    // Append tests to test suite
    append(static_cast<pfunction>(&TestGCTAResponse::test_response), "Test response");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_aeff), "Test effective area");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_psf), "Test PSF");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_npsf), "Test integrated PSF");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_irf_diffuse), "Test diffuse IRF");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_npred_diffuse), "Test diffuse IRF integration");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set CTA observation test methods
 ***************************************************************************/
void TestGCTAObservation::set(void)
{
    // Set test name
    name("GCTAObservation");

    // Append tests to test suite
    append(static_cast<pfunction>(&TestGCTAObservation::test_unbinned_obs), "Test unbinned observations");
    append(static_cast<pfunction>(&TestGCTAObservation::test_binned_obs), "Test binned observation");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set CTA optimizer test methods
 ***************************************************************************/
void TestGCTAOptimize::set(void)
{
    // Set test name
    name("CTA optimizers");

    // Append tests to test suite
    append(static_cast<pfunction>(&TestGCTAOptimize::test_unbinned_optimizer), "Test unbinned optimizer");
    append(static_cast<pfunction>(&TestGCTAOptimize::test_binned_optimizer), "Test binned optimizer");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test CTA Aeff computation
 ***************************************************************************/
void TestGCTAResponse::test_response_aeff(void)
{
    // Load response
    GCTAResponse rsp;
    rsp.caldb(cta_caldb);
    rsp.load(cta_irf);

    // Sum over effective area for control
    GEnergy      eng;
    double sum = 0.0;
    double ref = 154124059000.00006; //!< Adjust to actual value
    for (int i = 0; i < 30; ++i) {
        eng.TeV(pow(10.0, -1.7 + 0.1*double(i)));
        double aeff = rsp.aeff(0.0, 0.0, 0.0, 0.0, eng.log10TeV());
        //std::cout << eng << " " << eng.log10TeV() << " " << aeff << std::endl;
        sum += aeff;
    }
    test_value(sum, ref, 0.1, "Effective area verification");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test CTA psf computation
 *
 * The Psf computation is tested by integrating numerically the Psf
 * function. Integration is done in a rather simplistic way, by stepping
 * radially away from the centre. The integration is done for a set of
 * energies from 0.1-10 TeV.
 ***************************************************************************/
void TestGCTAResponse::test_response_psf(void)
{
    // Load response
    GCTAResponse rsp;
    rsp.caldb(cta_caldb);
    rsp.load(cta_irf);

    // Integrate Psf
    GEnergy eng;
    for (double e = 0.1; e < 10.0; e *= 2.0) {
        eng.TeV(e);
        double r     = 0.0;
        double dr    = 0.001;
        int    steps = int(1.0/dr);
        double sum   = 0.0;
        for (int i = 0; i < steps; ++i) {
            r   += dr;
            sum += rsp.psf(r*deg2rad, 0.0, 0.0, 0.0, 0.0, eng.log10TeV()) *
                   twopi * std::sin(r*deg2rad) * dr*deg2rad;
        }
        test_value(sum, 1.0, 0.001, "PSF integration for "+eng.print());
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test CTA npsf computation
 ***************************************************************************/
void TestGCTAResponse::test_response_npsf(void)
{
    // Setup CTA response
    GCTAResponse rsp;
    rsp.caldb(cta_caldb);
    rsp.load(cta_irf);

    // Setup npsf computation
    GSkyDir      srcDir;
    GEnergy      srcEng;
    GTime        srcTime;
    GCTAPointing pnt;
    GCTARoi      roi;
    GCTAInstDir  instDir;
    instDir.radec_deg(0.0, 0.0);
    roi.centre(instDir);
    roi.radius(2.0);
    srcEng.TeV(0.1);

    // Test PSF centred on ROI
    srcDir.radec_deg(0.0, 0.0);
    double npsf = rsp.npsf(srcDir, srcEng.log10TeV(), srcTime, pnt, roi);
    test_value(npsf, 1.0, 1.0e-3, "PSF(0,0) integration");

    // Test PSF offset but inside ROI
    srcDir.radec_deg(1.0, 1.0);
    npsf = rsp.npsf(srcDir, srcEng.log10TeV(), srcTime, pnt, roi);
    test_value(npsf, 1.0, 1.0e-3, "PSF(1,1) integration");

    // Test PSF outside and overlapping ROI
    srcDir.radec_deg(0.0, 2.0);
    npsf = rsp.npsf(srcDir, srcEng.log10TeV(), srcTime, pnt, roi);
    test_value(npsf, 0.492373, 1.0e-3, "PSF(0,2) integration");

    // Test PSF outside ROI
    srcDir.radec_deg(2.0, 2.0);
    npsf = rsp.npsf(srcDir, srcEng.log10TeV(), srcTime, pnt, roi);
    test_value(npsf, 0.0, 1.0e-3, "PSF(2,2) integration");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test CTA IRF computation for diffuse source model
 *
 * Tests the IRF computation for the diffuse source model. This is done
 * by calling the GCTAObservation::model method which in turn calls the
 * GCTAResponse::irf_diffuse method. The test is done for a small counts
 * map to keep the test executing reasonably fast.
 ***************************************************************************/
void TestGCTAResponse::test_response_irf_diffuse(void)
{
    // Set reference value
    double ref = 13803.800313356;

    // Set parameters
    double src_ra  = 201.3651;
    double src_dec = -43.0191;
    int    nebins  = 5;

    // Setup pointing on Cen A
    GSkyDir skyDir;
    skyDir.radec_deg(src_ra, src_dec);
    GCTAPointing pnt;
    pnt.dir(skyDir);

    // Setup skymap (10 energy layers)
    GSkymap map("CAR", "CEL", src_ra, src_dec, 0.5, 0.5, 10, 10, nebins);

    // Setup time interval
    GGti  gti;
    GTime tstart(0.0);
    GTime tstop(1800.0);
    gti.append(tstart, tstop);

    // Setup energy boundaries
    GEbounds ebounds;
    GEnergy  emin;
    GEnergy  emax;
    emin.TeV(0.1);
    emax.TeV(100.0);
    ebounds.setlog(emin, emax, nebins);

    // Setup event cube centered on Cen A
    GCTAEventCube cube(map, ebounds, gti);

    // Setup dummy CTA observation
    GCTAObservation obs;
    obs.ontime(1800.0);
    obs.livetime(1600.0);
    obs.deadc(1600.0/1800.0);
    obs.response(cta_irf, cta_caldb);
    obs.events(&cube);
    obs.pointing(pnt);

    // Load model for IRF computation
    GModels models(cta_rsp_xml);

    // Reset sum
    double sum = 0.0;

    // Iterate over all bins in event cube
    for (int i = 0; i < obs.events()->size(); ++i) {

        // Get event pointer
        const GEventBin* bin = (*(static_cast<const GEventCube*>(obs.events())))[i];

        // Get model and add to sum
        double model = obs.model(models, *bin, NULL) * bin->size();
        sum += model;

    }

    // Test sum
    test_value(sum, ref, 1.0e-5, "Diffuse IRF computation");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test CTA Npred computation
 *
 * Tests the Npred computation for the diffuse source model. This is done
 * by loading the model from the XML file and by calling the
 * GCTAObservation::npred method which in turn calls the
 * GCTAResponse::npred_diffuse method. The test takes a few seconds.
 ***************************************************************************/
void TestGCTAResponse::test_response_npred_diffuse(void)
{
    // Set reference value
    double ref = 11212.26274;

    // Set parameters
    double src_ra  = 201.3651;
    double src_dec = -43.0191;
    double roi_rad =   4.0;

    // Setup ROI centred on Cen A with a radius of 4 deg
    GCTARoi     roi;
    GCTAInstDir instDir;
    instDir.radec_deg(src_ra, src_dec);
    roi.centre(instDir);
    roi.radius(roi_rad);

    // Setup pointing on Cen A
    GSkyDir skyDir;
    skyDir.radec_deg(src_ra, src_dec);
    GCTAPointing pnt;
    pnt.dir(skyDir);

    // Setup dummy event list
    GGti     gti;
    GEbounds ebounds;
    GTime    tstart(0.0);
    GTime    tstop(1800.0);
    GEnergy  emin;
    GEnergy  emax;
    emin.TeV(0.1);
    emax.TeV(100.0);
    gti.append(tstart, tstop);
    ebounds.append(emin, emax);
    GCTAEventList events;
    events.roi(roi);
    events.gti(gti);
    events.ebounds(ebounds);

    // Setup dummy CTA observation
    GCTAObservation obs;
    obs.ontime(1800.0);
    obs.livetime(1600.0);
    obs.deadc(1600.0/1800.0);
    obs.response(cta_irf, cta_caldb);
    obs.events(&events);
    obs.pointing(pnt);

    // Load models for Npred computation
    GModels models(cta_rsp_xml);

    // Perform Npred computation
    double npred = obs.npred(models, NULL);

    // Test Npred
    test_value(npred, ref, 1.0e-5, "Diffuse Npred computation");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test CTA response handling
 ***************************************************************************/
void TestGCTAResponse::test_response(void)
{
    // Test CTA response loading
    test_try("Test CTA response loading");
    try {
        // Load response
        GCTAResponse rsp;
        rsp.caldb(cta_caldb);
        rsp.load(cta_irf);
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
 ***************************************************************************/
void TestGCTAObservation::test_unbinned_obs(void)
{
    // Set filenames
    const std::string file1 = "test_cta_obs_unbinned.xml";

    // Declare observations
    GObservations   obs;
    GCTAObservation run;

    // Load unbinned CTA observation
    test_try("Load unbinned CTA observation");
    try {
        run.load_unbinned(cta_events);
        run.response(cta_irf,cta_caldb);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Add observation (twice) to data
    test_try("Load unbinned CTA observation");
    try {
        obs.append(run);
        obs.append(run);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Loop over all events using iterators
    int num = 0;
    for (GObservations::iterator event = obs.begin(); event != obs.end(); ++event) {
        num++;
    }
    test_value(num, 8794, 1.0e-20, "Test observation iterator");

    // Loop over all events using iterator
    num = 0;
    GCTAEventList *ptr = static_cast<GCTAEventList*>(const_cast<GEvents*>(run.events()));
    for (GCTAEventList::iterator event = ptr->begin(); event != ptr->end(); ++event) {
        num++;
    }
    test_value(num, 4397, 1.0e-20, "Test event iterator");

    // Test XML loading
    test_try("Test XML loading");
    try {
        obs = GObservations(cta_unbin_xml);
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
 * @brief Test binned observation handling
 ***************************************************************************/
void TestGCTAObservation::test_binned_obs(void)
{
    // Set filenames
    const std::string file1 = "test_cta_obs_binned.xml";

    // Declare observations
    GObservations   obs;
    GCTAObservation run;

    // Load binned CTA observation
    test_try("Load unbinned CTA observation");
    try {
        run.load_binned(cta_cntmap);
        run.response(cta_irf,cta_caldb);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test XML loading
    test_try("Test XML loading");
    try {
        obs = GObservations(cta_bin_xml);
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
 * @brief Test unbinned optimizer
 ***************************************************************************/
void TestGCTAOptimize::test_unbinned_optimizer(void)
{
    // Declare observations
    GObservations   obs;
    GCTAObservation run;

    // Load unbinned CTA observation
    test_try("Load unbinned CTA observation");
    try {
        run.load_unbinned(cta_events);
        run.response(cta_irf,cta_caldb);
        obs.append(run);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Load models from XML file
    obs.models(cta_model_xml);

    // Perform LM optimization
    double fit_results[] = {83.6331, 0,
                            22.0145, 0,
                            5.656246512e-16, 1.91458426e-17,
                            -2.484100472, -0.02573396361,
                            300000, 0,
                            1, 0,
                            2.993705325, 0.03572658413,
                            6.490832107e-05, 1.749021094e-06,
                            -1.833584022, -0.01512223495,
                            1000000, 0,
                            1, 0};
    test_try("Perform LM optimization");
    try {
        GOptimizerLM opt;
        opt.max_iter(100);
        obs.optimize(opt);
        test_try_success();
        for (int i = 0, j = 0; i < obs.models().size(); ++i) {
            GModel* model = obs.models()[i];
            for (int k = 0; k < model->size(); ++k) {
                GModelPar& par  = (*model)[k];
                std::string msg = "Verify optimization result for " + par.print();
                test_value(par.real_value(), fit_results[j++], 5.0e-5, msg);
                test_value(par.real_error(), fit_results[j++], 5.0e-5, msg);
            }
        }
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Exit test
    return;

}


/***********************************************************************//**
 * @brief Test binned optimizer
 ***************************************************************************/
void TestGCTAOptimize::test_binned_optimizer(void)
{
    // Declare observations
    GObservations   obs;
    GCTAObservation run;

    // Load binned CTA observation
    test_try("Load binned CTA observation");
    try {
        run.load_binned(cta_cntmap);
        run.response(cta_irf,cta_caldb);
        obs.append(run);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Load models from XML file
    obs.models(cta_model_xml);

    // Perform LM optimization
    double fit_results[] = {83.6331, 0,
                            22.0145, 0,
                            5.616410411e-16, 1.904730785e-17,
                            -2.481781246, -0.02580905077,
                            300000, 0,
                            1, 0,
                            2.933677595, 0.06639644824,
                            6.550723074e-05, 1.945714239e-06,
                            -1.833781187, -0.0161464076,
                            1000000, 0,
                            1, 0};
    test_try("Perform LM optimization");
    try {
        GOptimizerLM opt;
        opt.max_iter(100);
        obs.optimize(opt);
        test_try_success();
        for (int i = 0, j = 0; i < obs.models().size(); ++i) {
            GModel* model = obs.models()[i];
            for (int k = 0; k < model->size(); ++k) {
                GModelPar& par  = (*model)[k];
                std::string msg = "Verify optimization result for " + par.print();
                test_value(par.real_value(), fit_results[j++], 5.0e-5, msg);
                test_value(par.real_error(), fit_results[j++], 5.0e-5, msg);
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
    GTestSuites testsuites("CTA instrument specific class testing");

    // Check if data directory exists
    bool has_data = (access(datadir.c_str(), R_OK) == 0);
    if (has_data) {
        std::string caldb = "CALDB="+cta_caldb;
        putenv((char*)caldb.c_str());
    }

    // Initially assume that we pass all tests
    bool success = true;

    // Create test suites and append them to the container
    TestGCTAResponse    rsp;
    TestGCTAObservation obs;
    TestGCTAOptimize    opt;
    testsuites.append(rsp);
    if (has_data) {
        testsuites.append(obs);
        testsuites.append(opt);
    }

    // Run the testsuites
    success = testsuites.run();

    // Save test report
    testsuites.save("reports/GCTA.xml");

    // Return success status
    return (success ? 0 : 1);
}
