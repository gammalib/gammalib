/***************************************************************************
 *                       test_CTA.cpp - Test CTA classes                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
#include <cmath>
#include <unistd.h>
#include "GCTALib.hpp"
#include "GTools.hpp"
#include "GNodeArray.hpp"
#include "test_CTA.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Constants __________________________________________________________ */
const std::string datadir          = PACKAGE_SOURCE"/inst/cta/test/data";
const std::string cta_caldb        = PACKAGE_SOURCE"/inst/cta/caldb";
const std::string cta_irf          = "cta_dummy_irf";
const std::string cta_events       = datadir+"/crab_events.fits.gz";
const std::string cta_cntmap       = datadir+"/crab_cntmap.fits.gz";
const std::string cta_cube_xml     = datadir+"/obs_cube.xml";
const std::string cta_bin_xml      = datadir+"/obs_binned.xml";
const std::string cta_unbin_xml    = datadir+"/obs_unbinned.xml";
const std::string cta_model_xml    = datadir+"/crab.xml";
const std::string cta_rsp_xml      = datadir+"/rsp_models.xml";
const std::string cta_cube_bgd_xml = datadir+"/cta_model_cube_bgd.xml";
const std::string cta_irf_bgd_xml  = datadir+"/cta_model_irf_bgd.xml";
const std::string cta_caldb_king   = PACKAGE_SOURCE"/inst/cta/caldb/data/cta/e/bcf/IFAE20120510_50h_King";
const std::string cta_irf_king     = "irf_file.fits";
const std::string cta_edisp_rmf    = PACKAGE_SOURCE"/inst/cta/test/caldb/dc1/rmf.fits";
const std::string cta_modbck_fit   = datadir+"/bg_test.fits";
const std::string cta_point_table  = datadir+"/crab_pointing.fits.gz";


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
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_psf_king), "Test King profile PSF");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_npsf), "Test integrated PSF");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_edisp), "Test energy dispersion");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_edispRMF), "Test energy dispersion RMF computation");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_irf_diffuse), "Test diffuse IRF");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_npred_diffuse), "Test diffuse IRF integration");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_expcube), "Test exposure cube");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_psfcube), "Test PSF cube");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGCTAResponse* TestGCTAResponse::clone(void) const
{
    // Clone test suite
    return new TestGCTAResponse(*this);
}


/***********************************************************************//**
 * @brief Set CTA model test methods
 ***************************************************************************/
void TestGCTAModel::set(void)
{
    // Set test name
    name("Test CTA models");

    // Append tests to test suite
    append(static_cast<pfunction>(&TestGCTAModel::test_model_cube_bgd), "Test CTA cube background model");
    append(static_cast<pfunction>(&TestGCTAModel::test_model_irf_bgd), "Test CTA IRF background model");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGCTAModel* TestGCTAModel::clone(void) const
{
    // Clone test suite
    return new TestGCTAModel(*this);
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
    append(static_cast<pfunction>(&TestGCTAObservation::test_cube_obs), "Test cube-style observation");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGCTAObservation* TestGCTAObservation::clone(void) const
{
    // Clone test suite
    return new TestGCTAObservation(*this);
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
    append(static_cast<pfunction>(&TestGCTAOptimize::test_cube_optimizer), "Test cube-style optimizer");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGCTAOptimize* TestGCTAOptimize::clone(void) const
{
    // Clone test suite
    return new TestGCTAOptimize(*this);
}




/***********************************************************************//**
 * @brief Set CTA pointing test methods
 ***************************************************************************/
void TestGCTAPointing::set(void)
{
    // Set test name
    name("GCTAPointing");

    // Append tests to test suite
    append(static_cast<pfunction>(&TestGCTAPointing::test_load_table), 
           "Test load pointing from table");

    append(static_cast<pfunction>(&TestGCTAPointing::test_interpolate_altaz),
           "Test alt/az interpolation given a time");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGCTAPointing* TestGCTAPointing::clone(void) const
{
    // Clone test suite
    return new TestGCTAPointing(*this);
}


/***********************************************************************//**
 * @brief Test ability to load a CTA pointing table
 ***************************************************************************/
void TestGCTAPointing::test_load_table(void)
{
    GCTAPointing pnt;
    pnt.load(cta_point_table);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test interpolation of alt/az pointing dir as a function of time
 ***************************************************************************/
void TestGCTAPointing::test_interpolate_altaz(void)
{

    GCTAObservation run;

    GCTAPointing pnt;
    pnt.load(cta_point_table);

    // Test an out-of bounds time
    test_try("Test an out-of bounds time");
    try {
        GTime time(155470378.7, "sec"); 
        GHorizDir dir = pnt.dir_horiz(time);
        test_try_failure();
    }
    catch (GException::out_of_range &exc) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Get a time that is somewhere in the run:
    GTime time(-128526000.0, "sec"); 
    GHorizDir dir = pnt.dir_horiz( time );

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test CTA response handling
 ***************************************************************************/
void TestGCTAResponse::test_response(void)
{
    // Test CTA IRF response constructor
    test_try("Test CTA IRF response constructor");
    try {
        GCTAResponseIrf rsp;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test CTA IRF response loading
    test_try("Test CTA IRF response loading");
    try {
        // Load response
        GCTAResponseIrf rsp;
        rsp.caldb(GCaldb(cta_caldb));
        rsp.load(cta_irf);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test CTA cube response constructors
    test_try("Test CTA cube response void constructor");
    try {
        GCTAResponseCube rsp;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    test_try("Test CTA cube response constructor");
    try {
        GCTAExposure exposure;
        GCTAMeanPsf psf;
        GCTAResponseCube rsp(exposure, psf);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test CTA Aeff computation
 ***************************************************************************/
void TestGCTAResponse::test_response_aeff(void)
{
    // Load response
    GCTAResponseIrf rsp;
    rsp.caldb(GCaldb(cta_caldb));
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
    GCTAResponseIrf rsp;
    rsp.caldb(GCaldb(cta_caldb));
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
            sum += rsp.psf(r * gammalib::deg2rad, 0.0, 0.0, 0.0, 0.0, eng.log10TeV()) *
                   gammalib::twopi * std::sin(r * gammalib::deg2rad) * dr *
                   gammalib::deg2rad;
        }
        test_value(sum, 1.0, 0.001, "PSF integration for "+eng.print());
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test CTA King profile psf computation
 *
 * The Psf computation is tested by integrating numerically the Psf
 * function. Integration is done in a rather simplistic way, by stepping
 * radially away from the centre. The integration is done for a set of
 * energies from 0.1-10 TeV.
 ***************************************************************************/
void TestGCTAResponse::test_response_psf_king(void)
{
    // Load response
    GCTAResponseIrf rsp;
    rsp.caldb(GCaldb(cta_caldb_king));
    rsp.load(cta_irf_king);

    // Integrate Psf
    GEnergy eng;
    for (double e = 0.1; e < 10.0; e *= 2.0) {
        eng.TeV(e);
        double r_max = rsp.psf()->delta_max(eng.log10TeV()) * gammalib::rad2deg;
        double r     = 0.0;
        double dr    = 0.0001;
        int    steps = int(r_max / dr);
        double sum   = 0.0;
        for (int i = 0; i < steps; ++i) {
            r   += dr;
            sum += rsp.psf(r * gammalib::deg2rad, 0.0, 0.0, 0.0, 0.0, eng.log10TeV()) *
                   gammalib::twopi * std::sin(r * gammalib::deg2rad) * dr *
                   gammalib::deg2rad;
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
    GCTAResponseIrf rsp;
    rsp.caldb(GCaldb(cta_caldb));
    rsp.load(cta_irf);

    // Setup npsf computation
    GSkyDir      srcDir;
    GEnergy      srcEng;
    GTime        srcTime;
    GCTAPointing pnt;
    GCTARoi      roi;
    GCTAInstDir  instDir;
    instDir.dir().radec_deg(0.0, 0.0);
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
 * @brief Test CTA Energy Dispersion computation
 *
 * The Energy Dispersion computation is tested by integrating numerically
 * the edisp function. Integration is done in a rather simplistic way, by
 * stepping through the energy range. The integration is done for a set of
 * true energies from 0.1-10 TeV.
 *
 * Note: test_energy_integration is separated out so that it may be more easily
 * called in other tests where it is useful (e.g.test_response_edispRMF and
 * future additions)
 ***************************************************************************/
void TestGCTAResponse::test_response_edisp(void)
{
    // Load response
	GCTAResponseIrf rsp;

	test_energy_integration(rsp);

	// Return
	return;
}


/***********************************************************************//**
 * @brief Test CTA Edisp RMF computation
 ***************************************************************************/
void TestGCTAResponse::test_response_edispRMF(void)
{
    // Load response
    GCTAResponseIrf rsp;
    rsp.caldb(GCaldb(cta_caldb));

    // Test Energy Dispersion
    test_energy_integration(rsp);

    // Test void constructor
    test_try("GRmf void constructor");
    try {
        GRmf rmf;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test source energy boundary constructor
    test_try("GRmf energy boundary constructor");
    try {
        GEbounds ebounds_src(10, GEnergy(0.1, "TeV"), GEnergy(10.0, "TeV"));
        GRmf     rmf(ebounds_src, ebounds_src);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test observed energy boundary constructor
    test_try("GRmf energy boundary constructor");
    try {
        GEbounds ebounds_obs(10, GEnergy(0.1, "TeV"), GEnergy(10.0, "TeV"));
        GRmf     rmf(ebounds_obs, ebounds_obs);
        test_try_success();
        }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    //
    GCTAEdispRmf edisp(cta_edisp_rmf);

    // Test if non-diagonal element (below diagonal) is zero
    test_value(edisp(std::log10(30.0),std::log10(1.0)), 0.0);

    // Test that diagonal element is non-zero
    //test_value(edisp(std::log10(30.0),std::log10(30.0)), 0.354952753 ,1e-9);

    // Test if non-diagonal element (above diagonal) is zero
    test_value(edisp(std::log10(1.0),std::log10(30.0)), 0.0);

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
    //double ref = 13803.800313356;
    //const double ref = 13803.5186374;
    const double ref = 13803.6932774; // After GWcs::solidangle improvement

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
    GEnergy  emin(0.1, "TeV");
    GEnergy  emax(100.0, "TeV");
    GEbounds ebounds(nebins, emin, emax);

    // Setup event cube centreed on Cen A
    GCTAEventCube cube(map, ebounds, gti);

    // Setup dummy CTA observation
    GCTAObservation obs;
    obs.ontime(1800.0);
    obs.livetime(1600.0);
    obs.deadc(1600.0/1800.0);
    obs.response(cta_irf, GCaldb(cta_caldb));
    obs.events(cube);
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
    //double ref = 11212.26274; // npred_spec precision of 1e-6
    //const double ref = 11212.437464; // npred_spec precision of 1e-5
    //const double ref = 11212.4370702; // After GWcs::solidangle improvement
    //const double ref = 12644.3391902; // After correcting for deadtime bug
    const double ref = 12643.6159142; // npred_spec precision of 1e-6

    // Set parameters
    double src_ra  = 201.3651;
    double src_dec = -43.0191;
    double roi_rad =   4.0;

    // Setup ROI centred on Cen A with a radius of 4 deg
    GCTARoi     roi;
    GCTAInstDir instDir;
    instDir.dir().radec_deg(src_ra, src_dec);
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
    GEnergy  emin(0.1, "TeV");
    GEnergy  emax(100.0, "TeV");
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
    obs.response(cta_irf, GCaldb(cta_caldb));
    obs.events(events);
    obs.pointing(pnt);

    // Load models for Npred computation
    GModels models(cta_rsp_xml);

    // Perform Npred computation
    double npred = obs.npred(models, NULL);

    // Test Npred
    test_value(npred, ref*1600.0/1800.0, 1.0e-5, "Diffuse Npred computation");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test exposure cube handling
 ***************************************************************************/
void TestGCTAResponse::test_response_expcube(void)
{
    // Test exposure cube constructors
    test_try("CTA exposure cube void constructor");
    try {
        GCTAExposure cube;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    test_try("CTA exposure cube map constructor");
    try {
        GEbounds     ebounds(20, GEnergy(0.1, "TeV"), GEnergy(100.0, "TeV"));
        GCTAExposure cube("CAR", "CEL", 83.63, 22.01, 0.02, 0.02, 200, 200, ebounds);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test set method
    GCTAObservation obs_cta;
    obs_cta.load(cta_cntmap);
    obs_cta.response(cta_irf, GCaldb(cta_caldb));
    GEbounds     ebounds(20, GEnergy(0.1, "TeV"), GEnergy(100.0, "TeV"));
    GCTAExposure cube("CAR", "CEL", 83.63, 22.01, 0.02, 0.02, 200, 200, ebounds);
    cube.set(obs_cta);
    cube.save("test_cta_expcube_one.fits", true);

    // Test fill method
    GObservations obs;
    obs_cta.id("000001");
    obs.append(obs_cta);
    obs_cta.id("000002");
    obs.append(obs_cta);
    cube.fill(obs);
    cube.save("test_cta_expcube_two.fits", true);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test PSF cube handling
 ***************************************************************************/
void TestGCTAResponse::test_response_psfcube(void)
{
    // Test PSF cube constructors
    test_try("CTA PSF cube void constructor");
    try {
        GCTAMeanPsf cube;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    test_try("CTA PSF cube map constructor");
    try {
        GEbounds    ebounds(20, GEnergy(0.1, "TeV"), GEnergy(100.0, "TeV"));
        GCTAMeanPsf cube("CAR", "CEL", 83.63, 22.01, 0.4, 0.4, 10, 10, ebounds, 0.1, 20);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test set method
    GCTAObservation obs_cta;
    obs_cta.load(cta_cntmap);
    obs_cta.response(cta_irf, GCaldb(cta_caldb));
    GEbounds    ebounds(20, GEnergy(0.1, "TeV"), GEnergy(100.0, "TeV"));
    GCTAMeanPsf cube("CAR", "CEL", 83.63, 22.01, 0.4, 0.4, 10, 10, ebounds, 0.1, 20);
    cube.set(obs_cta);
    cube.save("test_cta_psfcube_one.fits", true);

    // Test fill method
    GObservations obs;
    obs_cta.id("000001");
    obs.append(obs_cta);
    obs_cta.id("000002");
    obs.append(obs_cta);
    cube.fill(obs);
    cube.save("test_cta_psfcube_two.fits", true);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Utility function for energy dispersion tests
 ***************************************************************************/
void TestGCTAResponse::test_energy_integration(GCTAResponseIrf rsp) {

	// Load response
    rsp.caldb(GCaldb(cta_caldb));
    rsp.load(cta_irf);

	// Continue only if energy dispersion is available
	if (rsp.edisp() != NULL) {

	    // Loop over source energies (0.1 TeV -> 10.0 TeV)
	    for (double e_src = 0.1; e_src < 10.0; e_src *= 2.0) {

	        // Compute log10 of true energy
	        double log10_e_src = std::log10(e_src);

	        // Retrieve boundaries in observed energy
	        GEbounds ebounds  = rsp.edisp()->ebounds_obs(log10_e_src);
	        GEnergy  emin     = ebounds.emin();
	        GEnergy  emax     = ebounds.emax();
	        double   logE_min = std::log10(emin.TeV());
	        double   logE_max = std::log10(emax.TeV());

	        // Compute step size for numerical integration
	        const int steps = 1000;
	        double    dlogE = (logE_max-logE_min)/steps;

	        // Perform numerical integration by summing
	        double sum      = 0.0;
	        double logE_obs = logE_min;
	        for (int i = 0; i < steps; ++i) {
	            double dp_dlogE = (*rsp.edisp())(logE_obs, log10_e_src);
	            sum            += dp_dlogE * dlogE;
	            logE_obs       += dlogE;
	        }
	        GEnergy eng(e_src, "TeV");
	        test_value(sum, 1.0, 0.001, "Energy Dispersion integration for "+eng.print());

	        // And now in linear energies, test GCTAResponse::edisp() method
	        const int nE = 10000;
	        double E_min = emin.MeV();
	        double E_max = emax.MeV();
	        double dE    = (E_max-E_min)/nE;
	        double sum2  = 0.0;
	        double E_obs = E_min;
	        for (int i = 0; i < nE; ++i) {
	            double dp_dE = rsp.edisp(GEnergy(E_obs,"MeV"), 0.0, 0.0, 0.0, 0.0, log10_e_src);
	            sum2        += dp_dE * dE;
	            E_obs       += dE;
	        }
	        test_value(sum2, 1.0, 0.001, "Energy Dispersion integration for "+eng.print());
	    }

	} // endif: energy dispersion was available

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test CTA cube background
 ***************************************************************************/
void TestGCTAModel::test_model_cube_bgd(void)
{
    // Test CTA background constuctor
    test_try("Test CTA background constuctor");
    try {
        // Setup spectral model
        GEnergy energy(1.2, "TeV");
        GModelSpectralPlaw spectrum(1.0, -1.5, energy);

        // Load CTA event list
        GCTAObservation obs ;
        obs.load(cta_events);

        // Load background model
        GCTAModelCubeBackground bck(obs, cta_modbck_fit, spectrum);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Set reference value (~0.393469 result of a 2D Gaussian integrated 
    // within one sigma)
    const double ref = 1.0 - 1.0 / std::sqrt(std::exp(1.0));

    // Load models for Npred computation
    GModels models(cta_cube_bgd_xml);

    // Get the GCTAModelCubeBackground
    const GCTAModelCubeBackground* bck = dynamic_cast<const GCTAModelCubeBackground*>(models[0]);

    // Get the spectral and spatial components
    const GModelSpectralPlaw*       spec = dynamic_cast<const GModelSpectralPlaw*>(bck->spectral());
    const GModelSpatialRadialGauss* spat = dynamic_cast<const GModelSpatialRadialGauss*>(bck->spatial());

    // Get Integration centre for ROI position
    double src_ra  = spat->ra();
    double src_dec = spat->dec();
    double sigma   = spat->sigma();
    test_value(sigma, 1.0, 1e-7, "Input value from cta_modelbck.xml - file");

    // Set ROI to sigma of Gaussian
    double roi_rad = sigma;

    // Setup ROI centred on the Gaussian mean with a radius of 1sigma
    GCTARoi     roi;
    GCTAInstDir instDir;
    instDir.dir().radec_deg(src_ra, src_dec);
    roi.centre(instDir);
    roi.radius(roi_rad);

    // Setup pointing with the same centre as ROI
    GSkyDir skyDir;
    skyDir.radec_deg(src_ra, src_dec);
    GCTAPointing pnt;
    pnt.dir(skyDir);

    // Setup dummy event list
    GGti     gti;
    GEbounds ebounds;
    GTime    tstart(0.0);
    GTime    tstop(1800.0);
    GEnergy  emin(0.1, "TeV");
    GEnergy  emax(100.0, "TeV");
    gti.append(tstart, tstop);
    ebounds.append(emin, emax);
    GCTAEventList events;
    events.roi(roi);
    events.gti(gti);
    events.ebounds(ebounds);

    // Setup dummy CTA observation without deadtime
    GCTAObservation obs;
    obs.ontime(1800.0);
    obs.livetime(1800.0);
    obs.deadc(1.0);
    obs.response(cta_irf, GCaldb(cta_caldb));
    obs.events(events);
    obs.pointing(pnt);

    // Perform Npred computation
    double npred = bck->npred(spec->pivot(),tstart,obs);

    // Divide npred by spectral normalisation to true containment fraction
    npred /= spec->prefactor();

    // Test Npred against the reference value
    test_value(npred, ref , 1.0e-5,  "Npred computation for CTA background model");
	
	// Return
	return;
}


/***********************************************************************//**
 * @brief Test CTA IRF background model
 ***************************************************************************/
void TestGCTAModel::test_model_irf_bgd(void)
{
    // Test void constuctor
    test_try("Test void constuctor");
    try {
        GCTAModelIrfBackground model;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test XML constuctor
    test_try("Test XML constuctor");
    try {
        GXml xml(cta_irf_bgd_xml);
        const GXmlElement& lib = *xml.element("source_library", 0);
        const GXmlElement& src = *lib.element("source", 0);
        GCTAModelIrfBackground model(src);
        test_value(model["Prefactor"].value(), 1.0);
        test_value(model["Index"].value(), 0.0);
        test_value(model["PivotEnergy"].value(), 1.0e6);
        test_assert(model.is_constant(), "Model is expected to be constant.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test spectral constuctor
    test_try("Test spectral constuctor");
    try {
        GModelSpectralPlaw plaw(1.0, 0.0, GEnergy(1.0, "TeV"));
        GCTAModelIrfBackground model(plaw);
        test_value(model["Prefactor"].value(), 1.0);
        test_value(model["Index"].value(), 0.0);
        test_value(model["PivotEnergy"].value(), 1.0e6);
        test_assert(model.is_constant(), "Model is expected to be constant.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test XML loading of instrumental background
    GModels models(cta_irf_bgd_xml);
    GModel* model = models["My model"];
    test_value((*model)["Prefactor"].value(), 1.0);
    test_value((*model)["Index"].value(), 0.0);
    test_value((*model)["PivotEnergy"].value(), 1.0e6);
    test_assert(model->is_constant(), "Model is expected to be constant.");

    // Test XML saving and reloading of instrumental background
    models.save("test.xml");
    models.load("test.xml");
    model = models["My model"];
    test_value((*model)["Prefactor"].value(), 1.0);
    test_value((*model)["Index"].value(), 0.0);
    test_value((*model)["PivotEnergy"].value(), 1.0e6);
    test_assert(model->is_constant(), "Model is expected to be constant.");

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
        run.load(cta_events);
        run.response(cta_irf, GCaldb(cta_caldb));
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
    const GEvents* events = run.events();
    int num = 0;
    for (int i = 0; i < events->size(); ++i) {
        num++;
    }
    test_value(num, 4397, 1.0e-20, "Test event iterator");

    // Test XML loading and saving
    test_try("Test XML loading and saving");
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
    test_try("Load binned CTA observation");
    try {
        run.load(cta_cntmap);
        run.response(cta_irf, GCaldb(cta_caldb));
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test XML loading and saving
    test_try("Test XML loading and saving");
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
 * @brief Test cube-style observation handling
 ***************************************************************************/
void TestGCTAObservation::test_cube_obs(void)
{
    // Set filenames
    const std::string filename = "test_cta_obs_cube.xml";

    // Declare observations
    GObservations   obs;
    GCTAObservation run;

    // Load cube-style CTA observation
    test_try("Load cube-style CTA observation");
    try {
        run.load(cta_cntmap);
        //run.response(cta_irf, GCaldb(cta_caldb));
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test XML loading and saving
    test_try("Test XML loading and saving");
    try {
        obs = GObservations(cta_cube_xml);
        obs.save(filename);
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
        run.load(cta_events);
        run.response(cta_irf, GCaldb(cta_caldb));
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
                            -2.484100472, 0.02573396361,
                            300000, 0,
                            1, 0,
                            2.993705325, 0.03572658413,
                            6.490832107e-05, 1.749021094e-06,
                            -1.833584022, 0.01512223495,
                            1000000, 0,
                            1, 0};
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
        run.load(cta_cntmap);
        run.response(cta_irf, GCaldb(cta_caldb));
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
                            -2.48178, 0.02580905077,
                            300000, 0,
                            1, 0,
                            2.93368, 0.06639644824,
                            6.550723074e-05, 1.945714239e-06,
                            -1.833781187, 0.0160819,
                            1000000, 0,
                            1, 0};
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


/***********************************************************************//**
 * @brief Test cube optimizer
 ***************************************************************************/
void TestGCTAOptimize::test_cube_optimizer(void)
{
    // Load cube-style CTA observation
    GObservations obs(cta_cube_xml);

    // Load models from XML file
    obs.models(cta_model_xml);

    // Perform LM optimization
    double fit_results[] = {83.6331, 0,
                            22.0145, 0,
                            5.72212e-16, 2.01231e-17,
                            -2.48412, 0.0258234,
                            300000, 0,
                            1, 0,
                            2.93338, 0.0663668,
                            6.68923e-05, 1.96972e-06,
                            -1.83394, 0.0159896,
                            1000000, 0,
                            1, 0};
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
                test_value(par.value(), fit_results[j++], 1.0e-4, msg);
                test_value(par.error(), fit_results[j++], 1.0e-4, msg);
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
    // Allocate test suite container
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
    TestGCTAResponse        rsp;
    TestGCTAModel           model;
    TestGCTAOptimize        opt;
    TestGCTAObservation     obs;
    TestGCTAPointing        pnt;
    testsuites.append(rsp);
    if (has_data) {
    	testsuites.append(model);
        testsuites.append(opt);
        testsuites.append(obs);
        testsuites.append(pnt);
    }

    // Run the testsuites
    success = testsuites.run();

    // Save test report
    testsuites.save("reports/GCTA.xml");

    // Return success status
    return (success ? 0 : 1);
}
