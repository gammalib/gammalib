/***************************************************************************
 *                      test_CTA.cpp  -  test CTA classes                  *
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
 * @brief Testing of CTA classes
 * @author J. Knoedlseder
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
 * @brief Test CTA Aeff computation.
 ***************************************************************************/
void test_response_aeff(void)
{
    try {
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
        //printf("%30.5f\n", sum);
        if (std::abs(sum - ref) > 0.1) {
            std::cout << std::endl
                      << "TEST ERROR: Effective area differs from expected value"
                      << " (difference=" << std::abs(sum - ref) << ")"
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

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test CTA psf computation.
 ***************************************************************************/
void test_response_psf(void)
{
    // Test CTA Psf response
    try {
        // Load response
        GCTAResponse rsp;
        rsp.caldb(cta_caldb);
        rsp.load(cta_irf);

        // Integrate Psf
        GEnergy eng;
        for (double e = 0.1; e < 10.0; e*=2) {
            eng.TeV(e);
            double r     = 0.0;
            double dr    = 0.001;
            int    steps = int(1.0/dr);
            double sum   = 0.0;
            for (int i = 0; i < steps; ++i) {
                r   += dr;
                //obsDir.radec(0.0, r);
                sum += rsp.psf(r, 0.0, 0.0, 0.0, 0.0, eng.log10TeV()) *
                       twopi * std::sin(r*deg2rad) * dr;
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

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test CTA npsf computation.
 ***************************************************************************/
void test_response_npsf(void)
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

    // Try block to catch any problems in the computation
    try {

        // Test PSF centred on ROI
        srcDir.radec_deg(0.0, 0.0);
        double npsf = rsp.npsf(srcDir, srcEng.log10TeV(), srcTime, pnt, roi);
        double ref  = 1.0;
        if (std::abs(npsf - ref) > 1.0e-3) {
            std::cout << std::endl
                      << "TEST ERROR: Uncertainty in PSF(0,0) integration >0.1%"
                      << " (difference=" << std::abs(npsf - 1.0) << ")"
                      << std::endl;
            throw;
        }
        std::cout << ".";

        // Test PSF offset but inside ROI
        srcDir.radec_deg(1.0, 1.0);
        npsf = rsp.npsf(srcDir, srcEng.log10TeV(), srcTime, pnt, roi);
        ref  = 1.0;
        if (std::abs(npsf - ref) > 1.0e-3) {
            std::cout << std::endl
                      << "TEST ERROR: Uncertainty in PS(1,1) integration >0.1%"
                      << " (difference=" << std::abs(npsf - ref) << ")"
                      << std::endl;
            throw;
        }
        std::cout << ".";

        // Test PSF outside and overlapping ROI
        srcDir.radec_deg(0.0, 2.0);
        npsf = rsp.npsf(srcDir, srcEng.log10TeV(), srcTime, pnt, roi);
        ref  = 0.492373;
        if (std::abs(npsf - ref) > 1.0e-3) {
            std::cout << std::endl
                      << "TEST ERROR: Uncertainty in PS(0,2.1) integration >0.1%"
                      << " (difference=" << std::abs(npsf - ref) << ")"
                      << std::endl;
            throw;
        }
        std::cout << ".";

        // Test PSF outside ROI
        srcDir.radec_deg(2.0, 2.0);
        npsf = rsp.npsf(srcDir, srcEng.log10TeV(), srcTime, pnt, roi);
        ref  = 0.0;
        if (std::abs(npsf - ref) > 1.0e-3) {
            std::cout << std::endl
                      << "TEST ERROR: Uncertainty in PS(2,2) integration >0.1%"
                      << " (difference=" << std::abs(npsf - ref) << ")"
                      << std::endl;
            throw;
        }
        std::cout << ".";
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to compute GCTAResponse::npsf."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

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
void test_response_irf_diffuse(void)
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
    GTime tstart;
    GTime tstop;
    tstart.met(0.0);
    tstop.met(1800.0);
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

    // Extract CTA response
    //GResponse* rsp = obs.response();

    // Load model for IRF computation
    GModels models(cta_rsp_xml);

    // Try block to catch any problems in the computation
    try {

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
        if (std::abs(sum - ref) > 1.0e-5) {
            std::cout << std::endl
                      << "TEST ERROR: Diffuse IRF computation value differs"
                      << " from expectation by "
                      << sum - ref
                      << " (result-reference)."
                      << std::endl;
            throw;
        }

    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to perform diffuse IRF computation."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

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
void test_response_npred_diffuse(void)
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
    GTime    tstart;
    GTime    tstop;
    GEnergy  emin;
    GEnergy  emax;
    tstart.met(0.0);
    tstop.met(1800.0);
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

    // Extract CTA response
    //GResponse* rsp = obs.response();

    // Load models for Npred computation
    GModels models(cta_rsp_xml);

    // Try block to catch any problems in the computation
    try {

        // Perform Npred computation
        double npred = obs.npred(models, NULL);

        // Test Npred
        if (std::abs(npred - ref) > 1.0e-5) {
            std::cout << std::endl
                      << "TEST ERROR: Diffuse Npred computation value differs"
                      << " from expectation by "
                      << npred - ref
                      << " (result-reference)."
                      << std::endl;
            throw;
        }

    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to compute GCTAResponse::npred."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test CTA response handling.
 ***************************************************************************/
void test_response(void)
{
    // Dump header
    std::cout << "Test CTA response: ";

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
    test_response_aeff();

    // Test CTA Psf response
    test_response_psf();

    // Test GCTAResponse::npsf
    test_response_npsf();

    // Test GCTAResponse::irf_diffuse
    test_response_irf_diffuse();
    
    // Test GCTAResponse::npred_diffuse
    test_response_npred_diffuse();

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
    // Set filenames
    const std::string file1 = "test_cta_obs_unbinned.xml";

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
        if (num != 8794) {
            std::cout << std::endl
                      << "TEST ERROR: Wrong number of iterations in GObservations::iterator."
                      << " (excepted 8794, found " << num << ")" << std::endl;
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
        GCTAEventList *ptr = static_cast<GCTAEventList*>(const_cast<GEvents*>(run.events()));
        for (GCTAEventList::iterator event = ptr->begin(); event != ptr->end(); ++event) {
            //std::cout << *event->energy() << std::endl;
            num++;
        }
        if (num != 4397) {
            std::cout << std::endl
                      << "TEST ERROR: Wrong number of iterations in GCTAEventList::iterator."
                      << " (excepted 4397, found " << num << ")" << std::endl;
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

    // Test XML loading
    try {
        // Load observations
        obs = GObservations(cta_unbin_xml);
        
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
    std::cout << "Test CTA unbinned optimizer: ";

    // Number of observations in data
    //int nobs = 1;

    // Declare observations
    GObservations   obs;
    GCTAObservation run;

    // Load unbinned CTA observation
    try {
        // Load data and response and set ROI, energy range and time range
        // for analysis
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

    // Load models from XML file
    obs.models(cta_model_xml);

    // Perform LM optimization
    GOptimizerLM opt;
    try {
        opt.max_iter(100);
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
    std::cout << std::endl << opt << std::endl;
    std::cout << obs.models() << std::endl;

    // Plot final test success
    std::cout << " ok." << std::endl;

    // Exit test
    return;

}


/***********************************************************************//**
 * @brief Test unbinned observation handling.
 ***************************************************************************/
void test_binned_obs(void)
{
    // Set filenames
    const std::string file1 = "test_cta_obs_binned.xml";

    // Write header
    std::cout << "Test CTA binned observation handling: ";

    // Declare observations
    GObservations   obs;
    GCTAObservation run;

    // Load binned CTA observation
    try {
        run.load_binned(cta_cntmap);
        run.response(cta_irf,cta_caldb);
        //std::cout << run << std::endl;
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to load CTA run."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test XML loading
    try {
        // Load observations
        obs = GObservations(cta_bin_xml);
        
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

    // Notify final test success
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
    std::cout << "Test CTA binned optimizer: ";

    // Number of observations in data
    //int nobs = 1;

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

    // Load models from XML file
    obs.models(cta_model_xml);

    // Perform LM optimization
    GOptimizerLM opt;
    try {
        opt.max_iter(100);
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
    std::cout << std::endl << opt << std::endl;
    std::cout << obs.models() << std::endl;

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

    // Check if data directory exists
    bool has_data = (access(datadir.c_str(), R_OK) == 0);

    // Execute the tests
    test_response();

    // Execute tests requiring data
    if (has_data) {

        // Set CALDB environment variable
        std::string caldb = "CALDB="+cta_caldb;
        putenv((char*)caldb.c_str());

        // Execute tests
        test_unbinned_obs();
        test_binned_obs();
        test_unbinned_optimizer();
        test_binned_optimizer();
    }
    else {
        std::cout << "Skipped several tests since no test data have been found." << std::endl;
    }

    // Return
    return 0;
}
