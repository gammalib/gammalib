/***************************************************************************
 *                       test_CTA.cpp - Test CTA classes                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2018 by Juergen Knoedlseder                         *
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
#include <unistd.h>
#include <iostream>
#include <cmath>
#include <cstdlib>     // getenv
#include "GCTALib.hpp"
#include "GTools.hpp"
#include "GNodeArray.hpp"
#include "test_CTA.hpp"
#include "GCTAEdisp2D.hpp"
#include "GCTAEdispPerfTable.hpp"
#include "GCTAResponseTable.hpp"
#include "GMath.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Constants __________________________________________________________ */
const std::string datadir               = std::getenv("TEST_CTA_DATA");
const std::string caldbdir              = datadir + "/../caldb";
const std::string cta_caldb             = datadir + "/../../caldb";
const std::string cta_irf               = "cta_dummy_irf";
const std::string cta_events            = datadir+"/crab_events.fits";
const std::string cta_events_gti        = datadir+"/crab_events_gti.fits[EVENTS2]";
const std::string cta_cntmap            = datadir+"/crab_cntmap.fits";
const std::string cta_bin_xml           = datadir+"/obs_binned.xml";
const std::string cta_unbin_xml         = datadir+"/obs_unbinned.xml";
const std::string cta_model_xml         = datadir+"/crab.xml";
const std::string cta_rsp_xml           = datadir+"/rsp_models.xml";
const std::string cta_bgd_rad_gauss_xml = datadir+"/cta_model_bgd_rad_gauss.xml";
const std::string cta_cube_bgd_xml      = datadir+"/cta_model_cube_bgd.xml";
const std::string cta_irf_bgd_xml       = datadir+"/cta_model_irf_bgd.xml";
const std::string cta_aeff_bgd_xml      = datadir+"/cta_model_aeff_bgd.xml";
const std::string cta_caldb_king        = cta_caldb+"/data/cta/e/bcf/IFAE20120510_50h_King";
const std::string cta_irf_king          = "irf_file.fits";
const std::string cta_psf_table         = caldbdir+"/psf_table.fits[PSF_2D_TABLE]";
const std::string cta_perf_table        = caldbdir+"/cta_dummy_irf.dat";
const std::string cta_edisp_rmf         = caldbdir+"/dc1/rmf.fits";
const std::string cta_edisp_2D          = caldbdir+"/edisp_matrix.fits";
const std::string cta_bgd_3D            = caldbdir+"/edisp_matrix.fits";
const std::string cta_modbck_fit        = datadir+"/bg_test.fits";
const std::string cta_point_table       = datadir+"/crab_pointing.fits";

/* __ Test files for stacked analysis (based on Prod2::South_0.5h) _______ */
const std::string cta_stacked_xml       = datadir+"/stacked_obs.xml";
const std::string cta_stacked_model     = datadir+"/stacked_model.xml";
const std::string cta_stacked_cntcube   = datadir+"/stacked_cntcube.fits";
const std::string cta_stacked_expcube   = datadir+"/stacked_expcube.fits";
const std::string cta_stacked_psfcube   = datadir+"/stacked_psfcube.fits";
const std::string cta_stacked_edispcube = datadir+"/stacked_edispcube.fits";
const std::string cta_stacked_bkgcube   = datadir+"/stacked_bkgcube.fits";

/* __ Test files for On/Off analysis _____________________________________ */
const std::string cta_onoff_obs   = datadir+"/onoff_obs.xml";
const std::string cta_onoff_model = datadir+"/onoff_model.xml";
const std::string cta_onoff_onreg = datadir+"/onoff_region_on.reg";

/* __ Debug definitions __________________________________________________ */
//#define G_COMPUTE_REF_RATE_EBIN    //!< Compute reference rate


/***********************************************************************//**
 * @brief Set miscellaneous CTA classes test methods
 ***************************************************************************/
void TestGCTA::set(void)
{
    // Set test name
    name("Miscellaneous CTA classes");

    // Append GCTAInstDir tests to test suite
    append(static_cast<pfunction>(&TestGCTA::test_instdir),
           "Test GCTAInstDir class");

    // Append GCTAPointing tests to test suite
    append(static_cast<pfunction>(&TestGCTA::test_pointing_load_table),
           "Test loading of table into GCTAPointing instance");
    append(static_cast<pfunction>(&TestGCTA::test_pointing_interpolate_altaz),
           "Test alt/az interpolation given a time");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set CTA response test methods
 ***************************************************************************/
void TestGCTAResponse::set(void)
{
    // Set test name
    name("GCTAResponse");

    // Append tests to test suite
    append(static_cast<pfunction>(&TestGCTAResponse::test_response),
           "Test response");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_aeff),
           "Test effective area");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_psf),
           "Test Gaussian point spread function");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_psf_king),
           "Test King point spread function");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_psf_table),
           "Test Table point spread function");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_npsf),
           "Test integrated point spread function");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_irf_diffuse),
           "Test diffuse IRF");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_npred_diffuse),
           "Test diffuse IRF integration");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_edisp),
           "Test GCTAEdisp class");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_edisp_PerfTable),
           "Test GCTAEdispPerfTable class");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_edisp_RMF),
           "Test GCTAEdispRmf class");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_edisp_2D),
           "Test GCTAEdisp2D class");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_bgd_PerfTable),
           "Test GCTABackgroundPerfTable class");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_bgd_3D),
           "Test GCTABackground3D class");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_expcube),
           "Test GCTACubeExposure class");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_psfcube),
           "Test GCTACubePsf class");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_edispcube),
           "Test GCTACubeEdisp class");
    append(static_cast<pfunction>(&TestGCTAResponse::test_response_bkgcube),
           "Test GCTACubeBackground class");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set CTA model test methods
 ***************************************************************************/
void TestGCTAModel::set(void)
{
    // Set test name
    name("Test CTA models");

    // Append tests to test suite
    append(static_cast<pfunction>(&TestGCTAModel::test_model_bgd),
           "Test CTA background model");
    append(static_cast<pfunction>(&TestGCTAModel::test_model_cube_bgd),
           "Test CTA cube background model");
    append(static_cast<pfunction>(&TestGCTAModel::test_model_irf_bgd),
           "Test CTA IRF background model");
    append(static_cast<pfunction>(&TestGCTAModel::test_model_aeff_bgd),
           "Test CTA Aeff background model");

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
    append(static_cast<pfunction>(&TestGCTAObservation::test_event_bin),
           "Test event bin");
    append(static_cast<pfunction>(&TestGCTAObservation::test_event_cube),
           "Test event cube");
    append(static_cast<pfunction>(&TestGCTAObservation::test_unbinned_obs),
           "Test unbinned observations");
    append(static_cast<pfunction>(&TestGCTAObservation::test_binned_obs),
           "Test binned observation");
    append(static_cast<pfunction>(&TestGCTAObservation::test_stacked_obs),
           "Test stacked observation");
    append(static_cast<pfunction>(&TestGCTAObservation::test_onoff_obs),
           "Test On/Off observation");

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
    append(static_cast<pfunction>(&TestGCTAOptimize::test_unbinned_optimizer),
           "Test unbinned optimizer");
    append(static_cast<pfunction>(&TestGCTAOptimize::test_binned_optimizer),
           "Test binned optimizer");
    append(static_cast<pfunction>(&TestGCTAOptimize::test_stacked_optimizer),
           "Test stacked optimizer");
    append(static_cast<pfunction>(&TestGCTAOptimize::test_onoff_optimizer_cstat),
           "Test On/Off optimizer using CSTAT statistic");
    append(static_cast<pfunction>(&TestGCTAOptimize::test_onoff_optimizer_wstat),
           "Test On/Off optimizer using WSTAT statistic");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GCTAInstDir class
 ***************************************************************************/
void TestGCTA::test_instdir(void)
{
    // Check content of empty instance
    GCTAInstDir instdir1;

    // Test dir() method
    test_assert(!instdir1.has_dir(), "Test has_dir() method for empty instance");
    test_try("Test dir() method exception");
    try {
        GSkyDir dir = instdir1.dir();
        test_try_failure();
    }
    catch (GException::runtime_error &exc) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test detx() method
    test_assert(!instdir1.has_detx(), "Test has_detx() method for empty instance");
    test_try("Test has_detx() method exception");
    try {
        double detx = instdir1.detx();
        test_try_failure();
    }
    catch (GException::runtime_error &exc) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test dety() method
    test_assert(!instdir1.has_dety(), "Test has_dety() method for empty instance");
    test_try("Test has_dety() method exception");
    try {
        double dety = instdir1.dety();
        test_try_failure();
    }
    catch (GException::runtime_error &exc) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test theta() method
    test_try("Test theta() method exception");
    try {
        double dety = instdir1.theta();
        test_try_failure();
    }
    catch (GException::runtime_error &exc) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test phi() method
    test_try("Test phi() method exception");
    try {
        double dety = instdir1.phi();
        test_try_failure();
    }
    catch (GException::runtime_error &exc) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Set content
    double ra   = 83.6331 * gammalib::deg2rad;
    double dec  = 22.0145 * gammalib::deg2rad;
    double detx =     0.5 * gammalib::deg2rad;
    double dety =     1.5 * gammalib::deg2rad;
    GSkyDir dir;
    dir.radec(ra, dec);
    instdir1.dir(dir);
    instdir1.detx(detx);
    instdir1.dety(dety);

    // Check content of filled instance
    test_assert(instdir1.has_dir(), "Test has_dir() method for filled instance");
    test_assert(instdir1.has_detx(), "Test has_detx() method for filled instance");
    test_assert(instdir1.has_dety(), "Test has_dety() method for filled instance");
    test_value(instdir1.dir().ra(), ra, "Right Ascension of filled instance");
    test_value(instdir1.dir().dec(), dec, "Declination of filled instance");
    test_value(instdir1.detx(), detx, "DETX of filled instance");
    test_value(instdir1.dety(), dety, "DETY of filled instance");
    test_value(instdir1.theta(), 0.02759607852, "Theta of filled instance");
    test_value(instdir1.phi(), 1.24904577, "Phi of filled instance");

    // Check copying of instance
    GCTAInstDir instdir2(instdir1);
    test_assert(instdir2.has_dir(), "Test has_dir() method for copied instance");
    test_assert(instdir2.has_detx(), "Test has_detx() method for copied instance");
    test_assert(instdir2.has_dety(), "Test has_dety() method for copied instance");
    test_value(instdir2.dir().ra(), ra, "Right Ascension of copied instance");
    test_value(instdir2.dir().dec(), dec, "Declination of copied instance");
    test_value(instdir2.detx(), detx, "DETX of copied instance");
    test_value(instdir2.dety(), dety, "DETY of copied instance");
    test_value(instdir2.theta(), 0.02759607852, "Theta of copied instance");
    test_value(instdir2.phi(), 1.24904577, "Phi of copied instance");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test ability to load a CTA pointing table
 ***************************************************************************/
void TestGCTA::test_pointing_load_table(void)
{
    // Allocate classes
    GCTAPointing pnt;

    // Load pointing table
    pnt.load(cta_point_table);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test interpolation of alt/az pointing dir as a function of time
 ***************************************************************************/
void TestGCTA::test_pointing_interpolate_altaz(void)
{
    // Allocate classes
    GCTAObservation run;
    GCTAPointing    pnt;

    // Load pointing table
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
        GCTACubeExposure exposure;
        GCTACubePsf      psf;
        GCTACubeBackground      background;
        GCTAResponseCube rsp(exposure, psf, background);
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
    double  sum_check = 0.1;
    for (double e = 0.1; e < 10.0; e *= 2.0) {
        eng.TeV(e);
        double r     = 0.0;
        double dr    = 0.00001;
        int    steps = int(1.0/dr);
        double sum   = 0.0;
        double rcont = 0.0;
        for (int i = 0; i < steps; ++i) {
            r   += dr;
            sum += rsp.psf(r * gammalib::deg2rad, 0.0, 0.0, 0.0, 0.0, eng.log10TeV()) *
                   gammalib::twopi * std::sin(r * gammalib::deg2rad) * dr *
                   gammalib::deg2rad;
            
            // since 'sum' already totals to 1.0, its also a 'fraction',
            // which we can plug back into containment_radius(), to compare
            // with the origial radius 'r'
            if (sum > sum_check && sum < 1.0) {
                rcont      = rsp.psf()->containment_radius(sum, eng.log10TeV());
                sum_check += 0.1;
                std::string msg = "PSF containment radius for "+eng.print()+
                                  " and "+gammalib::str(sum*100.0)+
                                  "% containment";
                test_value(rcont, r*gammalib::deg2rad, 0.1, msg);
            }
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
    double  sum_check = 0.1;
    for (double e = 0.1; e < 10.0; e *= 2.0) {
        eng.TeV(e);
        double r_max = rsp.psf()->delta_max(eng.log10TeV()) * gammalib::rad2deg;
        double r     = 0.0;
        double dr    = 0.00001;
        int    steps = int(r_max / dr);
        double sum   = 0.0;
        double rcont = 0.0 ;
        for (int i = 0; i < steps; ++i) {
            r   += dr;
            sum += rsp.psf(r * gammalib::deg2rad, 0.0, 0.0, 0.0, 0.0, eng.log10TeV()) *
                   gammalib::twopi * std::sin(r * gammalib::deg2rad) * dr *
                   gammalib::deg2rad;
            
            // since 'sum' already totals to 1.0, its also a 'fraction',
            // which we can plug back into containment_radius(), to compare
            // with the origial radius 'r'
            if (sum > sum_check && sum < 1.0) {
                rcont      = rsp.psf()->containment_radius(sum, eng.log10TeV());
                sum_check += 0.1;
                std::string msg = "King PSF containment radius for "+eng.print()+
                                  " and "+gammalib::str(sum*100.0)+
                                  "% containment";
                test_value(rcont, r*gammalib::deg2rad, 0.1, msg);
            }
        }
        test_value(sum, 1.0, 0.005, "PSF integration for "+eng.print());
        
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test table point spread function
 *
 * The Psf computation is tested by integrating numerically the Psf
 * function. Integration is done in a rather simplistic way, by stepping
 * radially away from the centre. The integration is done for a set of
 * energies from 0.1-10 TeV.
 ***************************************************************************/
void TestGCTAResponse::test_response_psf_table(void)
{
    // Load point spread function
    GCTAPsfTable psf(cta_psf_table);

    // Integrate PSF
    GEnergy eng;
    double  sum_check = 0.1;
    for (double e = 0.1; e < 10.0; e *= 2.0) {
        eng.TeV(e);
        double r_max = psf.delta_max(eng.log10TeV()) * gammalib::rad2deg;
        double r     = 0.0;
        double dr    = 0.00001;
        int    steps = int(r_max / dr);
        double sum   = 0.0;
        double rcont = 0.0 ;
        for (int i = 0; i < steps; ++i) {
            r   += dr;
            sum += psf(r * gammalib::deg2rad, eng.log10TeV()) *
                   gammalib::twopi * std::sin(r * gammalib::deg2rad) * dr *
                   gammalib::deg2rad;

            // since 'sum' already totals to 1.0, its also a 'fraction',
            // which we can plug back into containment_radius(), to compare
            // with the origial radius 'r'
            if (sum > sum_check && sum < 1.0) {
                rcont      = psf.containment_radius(sum, eng.log10TeV());
                sum_check += 0.1;
                std::string msg = "Table PSF containment radius for "+eng.print()+
                                  " and "+gammalib::str(sum*100.0)+
                                  "% containment";
                test_value(rcont, r*gammalib::deg2rad, 0.1, msg);
            }
        }
        test_value(sum, 1.0, 0.01, "PSF integration for "+eng.print());

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
    GCTAInstDir  instDir(srcDir);
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
    //double ref = 13803.800313356;
    //const double ref = 13803.5186374;
    //const double ref = 13803.6932774; // after GWcs::solidangle improvement
    //const double ref = 13803.3453994; // after Irf computation improvements
    //const double ref = 13803.3461768; // using skymap flux() method
    //const double ref = 13803.5251338; // after correcting contains() method
    const double ref = 13803.058889257; // after increasing integration accuracy

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
    GSkyMap map("CAR", "CEL", src_ra, src_dec, 0.5, 0.5, 10, 10, nebins);

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
    test_value(sum, ref, 0.2, "Diffuse IRF computation");

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
    //const double ref = 11212.26274;   // npred_spec precision of 1e-6
    //const double ref = 11212.437464;  // npred_spec precision of 1e-5
    //const double ref = 11212.4370702; // After GWcs::solidangle improvement
    //const double ref = 12644.3391902; // After correcting for deadtime bug
    //const double ref = 12643.6159142; // npred_spec precision of 1e-6
    //const double ref = 10997.0851335; // After response computation change
    //const double ref = 11212.7161901; // npred_diffuse precision to iter=8
    //const double ref = 11238.6928076; // npred_diffuse precision to iter=9
    //const double ref = 11238.6933328; // using skymap flux() method
    const double ref = 11238.8390288; // after correcting contains() method

    // Set parameters
    double src_ra  = 201.3651;
    double src_dec = -43.0191;
    double roi_rad =   4.0;

    // Setup pointing on Cen A
    GSkyDir skyDir;
    skyDir.radec_deg(src_ra, src_dec);
    GCTAPointing pnt;
    pnt.dir(skyDir);

    // Setup ROI centred on Cen A with a radius of 4 deg
    GCTARoi     roi;
    GCTAInstDir instDir(skyDir);
    roi.centre(instDir);
    roi.radius(roi_rad);

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
    test_value(npred, ref, 0.2, "Diffuse Npred computation");

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
 ***************************************************************************/
void TestGCTAResponse::test_response_edisp(void)
{
    // Load response
    GCTAResponseIrf rsp;

	// Set performance table response
    rsp.caldb(GCaldb(cta_caldb));
    rsp.load(cta_irf);

    // Test energy dispersion integration
    test_response_edisp_integration(rsp);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test CTA performance table energy dispersion class
 ***************************************************************************/
void TestGCTAResponse::test_response_edisp_PerfTable(void)
{
    // Test GCTAEdispPerfTable void constructor
    GCTAEdispPerfTable edisp1;
    test_value(edisp1.filename(), "", "Test filename of void constructor");

    // Test GCTAEdispPerfTable filename constructor
    GCTAEdispPerfTable edisp2(cta_perf_table);
    test_value(edisp2.filename(), cta_perf_table,
               "Test filename of filename constructor");

    // Test load method
    edisp2.load(cta_perf_table);
    test_value(edisp2.filename(), cta_perf_table,
               "Test filename after loading perfromance table");

    // Test ebounds_obs method
    test_value(edisp2.ereco_bounds(GEnergy(0.1,"TeV")).emin().GeV(),  31.6795127349936);
    test_value(edisp2.ereco_bounds(GEnergy(0.1,"TeV")).emax().GeV(), 315.661420794325);

    // Test ebounds_src method
    test_value(edisp2.etrue_bounds(GEnergy(0.1,"TeV")).emin().GeV(),  31.6795127349936);
    test_value(edisp2.etrue_bounds(GEnergy(0.1,"TeV")).emax().GeV(), 315.661420794325);

    // Test if non-diagonal element (below diagonal) is zero
    test_value(edisp2(GEnergy(30.0,"TeV"), GEnergy(10.0,"TeV")), 0.0);

    // Test that diagonal element is non-zero
    test_value(edisp2(GEnergy(30.0,"TeV"), GEnergy(30.0,"TeV")), 2.12840930049511e-07);

    // Test if non-diagonal element (above diagonal) is zero
    test_value(edisp2(GEnergy(10.0,"TeV"), GEnergy(30.0,"TeV")), 0.0);


    // Test prob_erecobin method
    GEnergy etrue(1.0, "TeV");
    GEnergy emin(0.1, "TeV");
    GEnergy emax(10.0, "TeV");
    test_value(edisp2.prob_erecobin(emin, emax, etrue,  0.0), 1.0);
    test_value(edisp2.prob_erecobin(etrue, emax, etrue, 0.0), 0.5);
    test_value(edisp2.prob_erecobin(emin, etrue, etrue, 0.0), 0.5);

    // Test mc method
    GRan ran;
    test_value(edisp2.mc(ran, etrue).TeV(), 0.88667418590027);

    // Test normalisation
    test_edisp_integration(edisp2);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test CTA RMF energy dispersion class
 ***************************************************************************/
void TestGCTAResponse::test_response_edisp_RMF(void)
{
    // Test GCTAEdispRmf void constructor
    GCTAEdispRmf edisp1;
    test_value(edisp1.filename(), "", "Test filename of void constructor");

    // Test GCTAEdispRmf filename constructor
    GCTAEdispRmf edisp2(cta_edisp_rmf);
    test_value(edisp2.filename(), cta_edisp_rmf,
               "Test filename of filename constructor");

    // Test load method
    edisp2.load(cta_edisp_rmf);
    test_value(edisp2.filename(), cta_edisp_rmf,
               "Test filename after loading RMF file");

    // Test ebounds_obs method
    test_value(edisp2.ereco_bounds(GEnergy(10.0,"TeV")).emin().TeV(),  5.75439929962158);
    test_value(edisp2.ereco_bounds(GEnergy(10.0,"TeV")).emax().TeV(), 15.8489322662354);

    // Test ebounds_src method
    test_value(edisp2.etrue_bounds(GEnergy(10.0,"TeV")).emin().TeV(),  5.75439929962158);
    test_value(edisp2.etrue_bounds(GEnergy(10.0,"TeV")).emax().TeV(), 15.8489322662354);

    // Test if non-diagonal element (below diagonal) is zero
    test_value(edisp2(GEnergy(30.0,"TeV"), GEnergy(1.0,"TeV")), 0.0);

    // Test that diagonal element is non-zero
    test_value(edisp2(GEnergy(30.0,"TeV"), GEnergy(30.0,"TeV")), 1.07978131144085e-07);

    // Test if non-diagonal element (above diagonal) is zero
    test_value(edisp2(GEnergy(1.0,"TeV"), GEnergy(30.0,"TeV")), 0.0);

    // Test prob_erecobin method
    GEnergy etrue(1.0, "TeV");
    GEnergy emin(0.1, "TeV");
    GEnergy emax(10.0, "TeV");
    test_value(edisp2.prob_erecobin(emin, emax, etrue, 0.0),  1.00105708071354);
    test_value(edisp2.prob_erecobin(etrue, emax, etrue, 0.0), 0.50215659368534);
    test_value(edisp2.prob_erecobin(emin, etrue, etrue, 0.0), 0.498805373808917);

    // Test mc method
    GRan ran;
    test_value(edisp2.mc(ran, etrue).TeV(), 0.923744959981026);

    // Test normalisation
    test_edisp_integration(edisp2, 0.5, 10.0);

    // Return
    return;
}

/***********************************************************************//**
 * @brief Test CTA 2D energy dispersion class
 ***************************************************************************/
void TestGCTAResponse::test_response_edisp_2D(void)
{
    // Test GCTAEdisp2D void constructor
    GCTAEdisp2D edisp1;
    test_value(edisp1.filename(), "", "Test filename of void constructor");

    // Test GCTAEdisp2D filename constructor
    GCTAEdisp2D edisp2(cta_edisp_2D);
    test_value(edisp2.filename(), cta_edisp_2D,
               "Test filename of filename constructor");

    // Test load method
    edisp2.load(cta_edisp_2D);
    test_value(edisp2.filename(), cta_edisp_2D,
               "Test filename after loading response table file");

    // Test ebounds_obs method
    test_value(edisp2.ereco_bounds(GEnergy(10.0,"TeV")).emin().TeV(),  6.82816152379073);
    test_value(edisp2.ereco_bounds(GEnergy(10.0,"TeV")).emax().TeV(), 13.2011115155942);

    // Test ebounds_src method
    test_value(edisp2.etrue_bounds(GEnergy(10.0,"TeV")).emin().TeV(),  7.75156710957461);
    test_value(edisp2.etrue_bounds(GEnergy(10.0,"TeV")).emax().TeV(), 15.6834030551103);

    // Test if non-diagonal element (below diagonal) is zero
    test_value(edisp2(GEnergy(30.0,"TeV"), GEnergy(1.0,"TeV")), 0.0);

    // Test that diagonal element is non-zero
    test_value(edisp2(GEnergy(30.0,"TeV"), GEnergy(30.0,"TeV")), 1.74200329848547e-07);

    // Test if non-diagonal element (above diagonal) is zero
    test_value(edisp2(GEnergy(1.0,"TeV"), GEnergy(30.0,"TeV")), 0.0);

    // Test prob_erecobin method
    GEnergy etrue(1.0, "TeV");
    GEnergy emin(0.1, "TeV");
    GEnergy emax(10.0, "TeV");
    test_value(edisp2.prob_erecobin(emin, emax, etrue, 0.0),  1.00002101900895);
    test_value(edisp2.prob_erecobin(etrue, emax, etrue, 0.0), 0.480451962185927);
    test_value(edisp2.prob_erecobin(emin, etrue, etrue, 0.0), 0.519569056823021);

    // Test mc method
    GRan ran;
    test_value(edisp2.mc(ran, etrue).TeV(), 0.991300452406315);

    // Test normalisation
    test_edisp_integration(edisp2);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GCTABackgroundPerfTable class
 ***************************************************************************/
void TestGCTAResponse::test_response_bgd_PerfTable(void)
{
    // Set reference values
    const double ref_rate_below       = 0.000612295557409636;
    const double ref_rate_above       = 5.74072516267095e-11;
    const double ref_rate_onaxis      = 3.62782669235078e-05;
    const double ref_rate_offaxis     = 3.12350491143832e-05;
    const double ref_rate_ebin_within = 44.9181299290927;
    const double ref_rate_ebin_below  = 3676.68195289697;
    const double ref_rate_ebin_above  = 76.2328310391366;

    // Set energies for energy range for rate integration
    const GEnergy ebelow(10.0, "GeV");
    const GEnergy emin(1.0, "TeV");
    const GEnergy emax(10.0, "TeV");
    const GEnergy eabove(200.0, "TeV");

    // Test GCTABackgroundPerfTable constructor
    test_try("GCTABackgroundPerfTable constructor");
    try {
        GCTABackgroundPerfTable bgd1(cta_perf_table);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Set instrument direction
    GSkyDir     dir;
    GCTAInstDir instdir(dir, 0.0, 0.0);

    // Check content of empty instance
    GCTABackgroundPerfTable bgd1;
    test_value(bgd1.size(), 0, "Nodes of empty instance");
    test_value(bgd1.sigma(), 3.0, "Sigma value of empty instance");
    test_value(bgd1.filename(), "", "Filename of empty instance");
    test_value(bgd1(0.0, 0.0, 0.0), 0.0,
        "Onaxis operator value of empty instance");
    test_value(bgd1(0.0, 0.01, 0.02), 0.0,
        "Offaxis operator value of empty instance");
    test_value(bgd1.rate_ebin(instdir, emin, emax), 0.0,
        "Integrated rate value of empty instance");

    // Check content of filled instance
    GCTABackgroundPerfTable bgd2(cta_perf_table);
    test_value(bgd2.size(), 20, "Nodes of filled instance");
    test_value(bgd2.sigma(), 3.0, "Sigma value of filled instance");
    test_value(bgd2.filename(), cta_perf_table, "Filename of filled instance");
    test_value(bgd2(-2.0, 0.0, 0.0), ref_rate_below,
        "Operator value of filled instance extrapolated below lowest energy");
    test_value(bgd2(3.0, 0.0, 0.0), ref_rate_above,
        "Operator value of filled instance extrapolated above highest energy");
    test_value(bgd2(0.0, 0.0, 0.0), ref_rate_onaxis,
        "Onaxis operator value of filled instance");
    test_value(bgd2(0.0, 0.01, 0.02), ref_rate_offaxis,
        "Offaxis operator value of filled instance");
    test_value(bgd2.rate_ebin(instdir, emin, emax), ref_rate_ebin_within,
        "Integrated rate value of filled instance (within covered energies)");
    test_value(bgd2.rate_ebin(instdir, ebelow, emax), ref_rate_ebin_below,
        "Integrated rate value of filled instance (extrapolate below lowest energy)");
    test_value(bgd2.rate_ebin(instdir, emin, eabove), ref_rate_ebin_above,
        "Integrated rate value of filled instance (extrapolate above highest energy)");
    test_value(bgd2.rate_ebin(instdir, emin, emin), 0.0,
        "Integrated rate value of filled instance for zero interval");

    // Check content of copied instance
    GCTABackgroundPerfTable bgd3(bgd2);
    test_value(bgd3.size(), 20, "Nodes of copied instance");
    test_value(bgd3.sigma(), 3.0, "Sigma value of copied instance");
    test_value(bgd3.filename(), cta_perf_table, "Filename of copied instance");
    test_value(bgd3(-2.0, 0.0, 0.0), ref_rate_below,
        "Operator value of copied instance extrapolated below lowest energy");
    test_value(bgd3(3.0, 0.0, 0.0), ref_rate_above,
        "Operator value of copied instance extrapolated above highest energy");
    test_value(bgd3(0.0, 0.0, 0.0), ref_rate_onaxis,
        "Onaxis operator value of copied instance");
    test_value(bgd3(0.0, 0.01, 0.02), ref_rate_offaxis,
        "Offaxis operator value of copied instance");
    test_value(bgd3.rate_ebin(instdir, emin, emax), ref_rate_ebin_within,
        "Integrated rate value of copied instance");
    test_value(bgd3.rate_ebin(instdir, emin, emin), 0.0,
        "Integrated rate value of copied instance for zero interval");

    // Compute reference rate
    #if defined(G_COMPUTE_REF_RATE_EBIN)
    double  rate          = 0.0;
    int     nbins         = 10000;
    GEnergy E             = emin;
    GEnergy dE            = (eabove-emin)/nbins;
    for (int i = 0; i < nbins; ++i) {
        rate += bgd2(E.log10TeV(), 0.0, 0.0) * dE.MeV();
        E    += dE;
    }
    std::cout << rate << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GCTABackground3D class
 ***************************************************************************/
void TestGCTAResponse::test_response_bgd_3D(void)
{
    // Set reference values
    const double ref_rate_below       =     0.00402260284428378;
    const double ref_rate_above       =     0.0;
    const double ref_rate_onaxis      =     0.000133831125491202;
    const double ref_rate_offaxis     =     0.000128514485663256;
    const double ref_rate_ebin_within =   109.618252655233;
    const double ref_rate_ebin_below  = 14239.8196215465;
    const double ref_rate_ebin_above  =   114.826571413223;

    // Set energies for energy range for rate integration
    const GEnergy ebelow(10.0, "GeV");
    const GEnergy emin(1.0, "TeV");
    const GEnergy emax(10.0, "TeV");
    const GEnergy eabove(200.0, "TeV");

    // Set instrument direction
    GSkyDir     dir;
    GCTAInstDir instdir(dir, 0.0, 0.0);

    // Test GCTABackground3D constructor
    test_try("GCTABackground3D constructor");
    try {
        GCTABackground3D bgd1(cta_bgd_3D);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Check content of empty instance
    GCTABackground3D bgd1;
    test_value(bgd1.table().axes(), 0, "Response table axes of empty instance");
    test_value(bgd1.filename(), "", "Filename of empty instance");
    test_value(bgd1(0.0, 0.0, 0.0), 0.0,
        "Onaxis operator value of empty instance");
    test_value(bgd1(0.0, 0.01, 0.02), 0.0,
        "Offaxis operator value of empty instance");
    test_value(bgd1.rate_ebin(instdir, emin, emax), 0.0,
        "Integrated rate value of empty instance");

    // Check content of filled instance
    GCTABackground3D bgd2(cta_bgd_3D);
    test_value(bgd2.table().axes(), 3, "Response table axes of filled instance");
    test_value(bgd2.filename(), cta_bgd_3D, "Filename of filled instance");
    test_value(bgd2(-2.0, 0.0, 0.0), ref_rate_below,
        "Operator value of filled instance extrapolated below lowest energy");
    test_value(bgd2(3.0, 0.0, 0.0), ref_rate_above,
        "Operator value of filled instance extrapolated above highest energy");
    test_value(bgd2(0.0, 0.0, 0.0), ref_rate_onaxis,
        "Onaxis operator value of filled instance");
    test_value(bgd2(0.0, 0.01, 0.02), ref_rate_offaxis,
        "Offaxis operator value of filled instance");
    test_value(bgd2.rate_ebin(instdir, emin, emax), ref_rate_ebin_within,
        "Integrated rate value of filled instance");
    test_value(bgd2.rate_ebin(instdir, ebelow, emax), ref_rate_ebin_below,
        "Integrated rate value of filled instance (extrapolate below lowest energy)");
    test_value(bgd2.rate_ebin(instdir, emin, eabove), ref_rate_ebin_above,
        "Integrated rate value of filled instance (extrapolate above highest energy)");
    test_value(bgd2.rate_ebin(instdir, emin, emin), 0.0,
        "Integrated rate value of filled instance for zero interval");

    // Check content of copied instance
    GCTABackground3D bgd3(bgd2);
    test_value(bgd3.table().axes(), 3, "Response table axes of copied instance");
    test_value(bgd3.filename(), cta_bgd_3D, "Filename of copied instance");
    test_value(bgd3(-2.0, 0.0, 0.0), ref_rate_below,
        "Operator value of copied instance extrapolated below lowest energy");
    test_value(bgd3(3.0, 0.0, 0.0), ref_rate_above,
        "Operator value of copied instance extrapolated above highest energy");
    test_value(bgd3(0.0, 0.0, 0.0), ref_rate_onaxis,
        "Onaxis operator value of copied instance");
    test_value(bgd3(0.0, 0.01, 0.02), ref_rate_offaxis,
        "Offaxis operator value of copied instance");
    test_value(bgd3.rate_ebin(instdir, emin, emax), ref_rate_ebin_within,
        "Integrated rate value of copied instance");
    test_value(bgd3.rate_ebin(instdir, emin, emin), 0.0,
        "Integrated rate value of copied instance for zero interval");

    // Compute reference rate
    #if defined(G_COMPUTE_REF_RATE_EBIN)
    double  rate          = 0.0;
    int     nbins         = 10000;
    GEnergy E             = emin;
    GEnergy dE            = (eabove-E)/nbins;
    for (int i = 0; i < nbins; ++i) {
        rate += bgd2(E.log10TeV(), 0.0, 0.0) * dE.MeV();
        E    += dE;
    }
    std::cout << rate << std::endl;
    #endif

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
        GCTACubeExposure cube;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    test_try("CTA exposure cube map constructor");
    try {
        GEnergies        energies(20, GEnergy(0.1, "TeV"), GEnergy(100.0, "TeV"));
        GCTACubeExposure cube("CAR", "CEL", 83.63, 22.01, 0.02, 0.02, 200, 200, energies);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test set method
    GCTAObservation obs_cta;
    obs_cta.load(cta_events);
    obs_cta.response(cta_irf, GCaldb(cta_caldb));
    GEnergies        energies(20, GEnergy(0.1, "TeV"), GEnergy(100.0, "TeV"));
    GCTACubeExposure cube("CAR", "CEL", 83.63, 22.01, 0.02, 0.02, 200, 200, energies);
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
        GCTACubePsf cube;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    test_try("CTA PSF cube map constructor");
    try {
        GEnergies   energies(20, GEnergy(0.1, "TeV"), GEnergy(100.0, "TeV"));
        GCTACubePsf cube("CAR", "CEL", 83.63, 22.01, 0.4, 0.4, 10, 10, energies, 0.1, 20);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test set method
    GCTAObservation obs_cta;
    obs_cta.load(cta_events);
    obs_cta.response(cta_irf, GCaldb(cta_caldb));
    GEnergies   energies(20, GEnergy(0.1, "TeV"), GEnergy(100.0, "TeV"));
    GCTACubePsf cube("CAR", "CEL", 83.63, 22.01, 0.4, 0.4, 10, 10, energies, 0.1, 20);
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
 * @brief Test background cube handling
 ***************************************************************************/
void TestGCTAResponse::test_response_bkgcube(void)
{
    // Test background cube constructors
    test_try("CTA background cube void constructor");
    try {
        GCTACubeBackground cube;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    GCTACubeBackground cube;
    cube.load(cta_stacked_bkgcube);
    cube.save("test_cta_bkgcube.fits", true);

    // Return
    return;
}

/***********************************************************************//**
 * @brief Test energy dispersion cube handling
 ***************************************************************************/
void TestGCTAResponse::test_response_edispcube(void)
{
    // Test energy dispersion cube constructors
    test_try("CTA energy dispersion cube void constructor");
    try {
    	GCTACubeEdisp cube;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    test_try("CTA energy dispersion cube map constructor");
    try {
        GEnergies     energies(20, GEnergy(0.1, "TeV"), GEnergy(100.0, "TeV"));
        GCTACubeEdisp cube("CAR", "CEL", 83.63, 22.01, 0.4, 0.4, 10, 10,
                           energies, 2.0, 20);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Load IRF for energy dispersion computation
    GCTAObservation obs_cta;
    obs_cta.load(cta_events);
    obs_cta.response(cta_irf, GCaldb(cta_caldb));

    // Initialise energy dispersion cube
    GEnergies     energies(100, GEnergy(0.1, "TeV"), GEnergy(100.0, "TeV"));
    GCTACubeEdisp cube("CAR", "CEL", 83.63, 22.01, 0.4, 0.4, 10, 10,
                       energies, 3.0, 100);

    // Test setting of observation
    cube.set(obs_cta);

    // Test energy dispersion integration
    test_edispcube_integration(cube, 0.1, 10.0);

    // Test filling of single observation
    GObservations obs;
    obs_cta.id("000001");
    obs.append(obs_cta);
    cube.fill(obs);

    // Test energy dispersion integration
    test_edispcube_integration(cube, 0.1, 10.0);

    // Test fill of two observations
    obs_cta.id("000002");
    obs.append(obs_cta);
    cube.fill(obs);

    // Test energy dispersion integration
    test_edispcube_integration(cube, 0.1, 10.0);

    // Test saving
    cube.save("test_cta_edispcube.fits", true);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test response cache handling
 ***************************************************************************/
void TestGCTAResponse::test_response_cache(void)
{
    // Set some energies
    const GEnergy eng1(1.0, "TeV");
    const GEnergy eng2(2.0, "TeV");
    const GEnergy eng3(3.0, "TeV");

    // Set some instrument directions
    GSkyDir skydir1;
    GSkyDir skydir2;
    skydir1.radec_deg(1.0,10.0);
    skydir2.radec_deg(2.0,20.0);
    const GCTAInstDir dir1(skydir1, 0.0, 0.0);
    const GCTAInstDir dir2(skydir2, 0.0, 0.0);

    // Initialise results
    double value = 0.0;
    bool   flag  = false;

    // Test empty response cache
    GCTAResponseCache cache1;
    test_assert(cache1.is_empty(), "Test is_empty() method for empty cache");
	test_value(cache1.size(), 0, "Test size() method for empty cache");
	test_value(cache1.ndirs(), 0, "Test ndirs() method for empty cache");
	test_value(cache1.nerecos(), 0, "Test nerecos() method for empty cache");
	test_value(cache1.netrues(), 0, "Test netrues() method for empty cache");
    flag = cache1.contains("Crab", eng1, eng1, &value);
    test_assert(!flag, "Test contains() method flag for empty cache");

    // Test filled response cube (one model, two ereco=etrue)
    GCTAResponseCache cache2;
    cache2.set("Crab", eng1, eng1, 1.0);
    cache2.set("Crab", eng2, eng2, 2.0);
    test_assert(!cache2.is_empty(), "Test is_empty() method for filled cache");
	test_value(cache2.size(), 2, "Test size() method for filled cache");
	test_value(cache2.ndirs(), 1, "Test ndirs() method for filled cache");
	test_value(cache2.nerecos(), 2, "Test nerecos() method for filled cache");
	test_value(cache2.netrues(), 2, "Test netrues() method for filled cache");
    flag = cache2.contains("Crab", eng1, eng1, &value);
    test_assert(flag, "Test contains() method flag for filled cache");
	test_value(value, 1.0, "Test contains() method value for filled cache");
    flag = cache2.contains("Crab", eng2, eng2, &value);
    test_assert(flag, "Test contains() method flag for filled cache");
	test_value(value, 2.0, "Test contains() method value for filled cache");
    flag = cache2.contains("Crab", eng2, eng1, &value);
    test_assert(!flag, "Test contains() method flag for non-existing true energy");
    flag = cache2.contains("Crab", eng1, eng2, &value);
    test_assert(!flag, "Test contains() method flag for non-existing reco. energy");
    flag = cache2.contains("Vela", eng2, eng2, &value);
    test_assert(!flag, "Test contains() method flag for non-existing name");

    // Test clearing of cache
    cache2.clear();
    test_assert(cache2.is_empty(), "Test is_empty() method for cleared cache");
	test_value(cache2.size(), 0, "Test size() method for cleared cache");
	test_value(cache2.ndirs(), 0, "Test ndirs() method for empty cache");
	test_value(cache2.nerecos(), 0, "Test nerecos() method for cleared cache");
	test_value(cache2.netrues(), 0, "Test netrues() method for cleared cache");
    flag = cache2.contains("Crab", eng1, eng1, &value);
    test_assert(!flag, "Test contains() method flag for cleared cache");

    // Test filled response cube (two models, two ereco, two etrue per ereco)
    GCTAResponseCache cache3;
    cache3.set("Crab", eng1, eng1, 1.0);
    cache3.set("Crab", eng1, eng2, 2.0);
    cache3.set("Crab", eng2, eng1, 3.0);
    cache3.set("Crab", eng2, eng2, 4.0);
    cache3.set("Vela", eng1, eng1, 10.0);
    cache3.set("Vela", eng1, eng2, 11.0);
    cache3.set("Vela", eng2, eng1, 12.0);
    cache3.set("Vela", eng2, eng2, 13.0);
    test_assert(!cache3.is_empty(), "Test is_empty() method for filled cache");
	test_value(cache3.size(), 8, "Test size() method for filled cache");
	test_value(cache3.ndirs(), 1, "Test ndirs() method for filled cache");
	test_value(cache3.nerecos(), 4, "Test nerecos() method for filled cache");
	test_value(cache3.netrues(), 8, "Test netrues() method for filled cache");

    // Test Crab cache
    flag = cache3.contains("Crab", eng1, eng1, &value);
    test_assert(flag, "Test contains() method flag for filled cache");
	test_value(value, 1.0, "Test contains() method value for filled cache");
    flag = cache3.contains("Crab", eng1, eng2, &value);
    test_assert(flag, "Test contains() method flag for filled cache");
	test_value(value, 2.0, "Test contains() method value for filled cache");
    flag = cache3.contains("Crab", eng2, eng1, &value);
    test_assert(flag, "Test contains() method flag for filled cache");
	test_value(value, 3.0, "Test contains() method value for filled cache");
    flag = cache3.contains("Crab", eng2, eng2, &value);
    test_assert(flag, "Test contains() method flag for filled cache");
	test_value(value, 4.0, "Test contains() method value for filled cache");
    flag = cache3.contains("Crab", eng3, eng1, &value);
    test_assert(!flag, "Test contains() method flag for non-existing true energy");
    flag = cache3.contains("Crab", eng1, eng3, &value);
    test_assert(!flag, "Test contains() method flag for non-existing reco. energy");

    // Test Vela cache
    flag = cache3.contains("Vela", eng1, eng1, &value);
    test_assert(flag, "Test contains() method flag for filled cache");
	test_value(value, 10.0, "Test contains() method value for filled cache");
    flag = cache3.contains("Vela", eng1, eng2, &value);
    test_assert(flag, "Test contains() method flag for filled cache");
	test_value(value, 11.0, "Test contains() method value for filled cache");
    flag = cache3.contains("Vela", eng2, eng1, &value);
    test_assert(flag, "Test contains() method flag for filled cache");
	test_value(value, 12.0, "Test contains() method value for filled cache");
    flag = cache3.contains("Vela", eng2, eng2, &value);
    test_assert(flag, "Test contains() method flag for filled cache");
	test_value(value, 13.0, "Test contains() method value for filled cache");
    flag = cache3.contains("Vela", eng3, eng1, &value);
    test_assert(!flag, "Test contains() method flag for non-existing true energy");
    flag = cache3.contains("Vela", eng1, eng3, &value);
    test_assert(!flag, "Test contains() method flag for non-existing reco. energy");

    // Test non-existing name
    flag = cache3.contains("Orion", eng1, eng1, &value);
    test_assert(!flag, "Test contains() method flag for non-existing name");

    // Test filled response cube (two models, two dir, two ereco, two etrue per ereco)
    GCTAResponseCache cache4;
    cache4.set("Crab", dir1, eng1, eng1, 1.0);
    cache4.set("Crab", dir1, eng1, eng2, 2.0);
    cache4.set("Crab", dir1, eng2, eng1, 3.0);
    cache4.set("Crab", dir1, eng2, eng2, 4.0);
    cache4.set("Crab", dir2, eng1, eng1, 5.0);
    cache4.set("Crab", dir2, eng1, eng2, 6.0);
    cache4.set("Crab", dir2, eng2, eng1, 7.0);
    cache4.set("Crab", dir2, eng2, eng2, 8.0);
    test_assert(!cache4.is_empty(), "Test is_empty() method for filled cache");
	test_value(cache4.size(), 8, "Test size() method for filled cache");
	test_value(cache4.ndirs(), 2, "Test ndirs() method for filled cache");
	test_value(cache4.nerecos(), 4, "Test nerecos() method for filled cache");
	test_value(cache4.netrues(), 8, "Test netrues() method for filled cache");

    // Test Crab cache for dir1
    flag = cache4.contains("Crab", dir1, eng1, eng1, &value);
    test_assert(flag, "Test contains() method flag for filled cache");
	test_value(value, 1.0, "Test contains() method value for filled cache");
    flag = cache4.contains("Crab", dir1, eng1, eng2, &value);
    test_assert(flag, "Test contains() method flag for filled cache");
	test_value(value, 2.0, "Test contains() method value for filled cache");
    flag = cache4.contains("Crab", dir1, eng2, eng1, &value);
    test_assert(flag, "Test contains() method flag for filled cache");
	test_value(value, 3.0, "Test contains() method value for filled cache");
    flag = cache4.contains("Crab", dir1, eng2, eng2, &value);
    test_assert(flag, "Test contains() method flag for filled cache");
	test_value(value, 4.0, "Test contains() method value for filled cache");
    flag = cache4.contains("Crab", dir1, eng3, eng1, &value);
    test_assert(!flag, "Test contains() method flag for non-existing true energy");
    flag = cache4.contains("Crab", dir1, eng1, eng3, &value);
    test_assert(!flag, "Test contains() method flag for non-existing reco. energy");

    // Test Crab cache for dir2
    flag = cache4.contains("Crab", dir2, eng1, eng1, &value);
    test_assert(flag, "Test contains() method flag for filled cache");
	test_value(value, 5.0, "Test contains() method value for filled cache");
    flag = cache4.contains("Crab", dir2, eng1, eng2, &value);
    test_assert(flag, "Test contains() method flag for filled cache");
	test_value(value, 6.0, "Test contains() method value for filled cache");
    flag = cache4.contains("Crab", dir2, eng2, eng1, &value);
    test_assert(flag, "Test contains() method flag for filled cache");
	test_value(value, 7.0, "Test contains() method value for filled cache");
    flag = cache4.contains("Crab", dir2, eng2, eng2, &value);
    test_assert(flag, "Test contains() method flag for filled cache");
	test_value(value, 8.0, "Test contains() method value for filled cache");
    flag = cache4.contains("Crab", dir2, eng3, eng1, &value);
    test_assert(!flag, "Test contains() method flag for non-existing true energy");
    flag = cache4.contains("Crab", dir2, eng1, eng3, &value);
    test_assert(!flag, "Test contains() method flag for non-existing reco. energy");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Utility function for energy dispersion tests
 *
 * @param[in] rsp Response.
 * @param[in] e_src_min Minimum true energy (TeV).
 * @param[in] e_src_max Maximum true energy (TeV).
 ***************************************************************************/
void TestGCTAResponse::test_response_edisp_integration(const GCTAResponseIrf& rsp,
                                                       const double&          e_src_min,
                                                       const double&          e_src_max)
{
	// Continue only if energy dispersion is available
	if (rsp.edisp() != NULL) {

	    // Loop over source energies
        for (double e_src = e_src_min; e_src < e_src_max; e_src *= 2.0) {

	        // Compute log10 of true energy
            GEnergy etrue(e_src,"TeV");

	        // Retrieve observed energy boundaries
	        GEbounds ebounds = rsp.edisp()->ereco_bounds(etrue);
	        double   emin    = ebounds.emin().MeV();
	        double   emax    = ebounds.emax().MeV();

	        // Perform numerical integration by summing
	        const int nE = 10000;
	        double dE    = (emax-emin)/nE;
	        double sum1  = 0.0;
	        double sum2  = 0.0;
	        double E_obs = emin;
	        for (int i = 0; i < nE; ++i, E_obs += dE) {
                GEnergy ereco(E_obs,"MeV");
	            double dp1  = (*rsp.edisp())(ereco, etrue);
	            double dp2  = rsp.edisp(ereco, etrue, 0.0, 0.0, 0.0, 0.0);
	            sum1       += dp1 * dE;
	            sum2       += dp2 * dE;
	        }
	        test_value(sum1, 1.0, 0.001, "Energy dispersion integration using "+
                                         rsp.edisp()->classname()+"::operator() "
                                         "for "+etrue.print());
	        test_value(sum2, 1.0, 0.001, "Energy dispersion integration using "
                                         "GCTAResponseIrf::edisp() for "+etrue.print());
	    }

	} // endif: energy dispersion was available

    // Return
    return;
}


/***********************************************************************//**
 * @brief Utility function for energy dispersion integration tests
 *
 * @param[in] edisp Energy dispersion instance.
 * @param[in] e_src_min Minimum true energy (TeV).
 * @param[in] e_src_max Maximum true energy (TeV).
 *
 * Tests whether the energy dispersion integrated over observed energies
 * gives unity.
 ***************************************************************************/
void TestGCTAResponse::test_edisp_integration(const GCTAEdisp&   edisp,
                                              const double&      e_src_min,
                                              const double&      e_src_max)
{
    // Loop over source energies
    for (double e_src = e_src_min; e_src <= e_src_max; e_src *= 2.0) {

        // Compute log10 of true energy
        GEnergy etrue(e_src,"TeV");

        // Retrieve observed energy boundaries
        GEbounds ebounds = edisp.ereco_bounds(etrue);

        // Skip this energy if the boundaries are empty
        if (ebounds.is_empty()) {
            continue;
        }

        // Perform numerical integration by summing
        const int nE = 10000;
        double emin  = ebounds.emin().MeV();
        double emax  = ebounds.emax().MeV();
        double dE    = (emax-emin)/nE;
        double sum   = 0.0;
        double E_obs = emin;
        for (int i = 0; i < nE; ++i, E_obs += dE) {
            GEnergy ereco(E_obs,"MeV");
            double dp  = edisp(ereco, etrue);
            sum       += dp * dE;
        }
        test_value(sum, 1.0, 0.001, "Energy dispersion integration using "+
                                    edisp.classname()+"::operator() for "+
                                    etrue.print());
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Utility function for energy dispersion cube integration tests
 *
 * @param[in] edisp Energy dispersion cube instance.
 * @param[in] e_src_min Minimum true energy (TeV).
 * @param[in] e_src_max Maximum true energy (TeV).
 *
 * Tests whether the energy dispersion cube integrated over observed
 * energies gives unity.
 ***************************************************************************/
void TestGCTAResponse::test_edispcube_integration(const GCTACubeEdisp& edisp,
                                                  const double&        e_src_min,
                                                  const double&        e_src_max)
{
    // Set sky direction
    GSkyDir dir;
    dir.radec_deg(83.63, 22.01);

    // Test energy dispersion cube integration
    for (double e_src = e_src_min; e_src <= e_src_max; e_src *= 2.0) {

        // Compute log10 of true energy
        GEnergy etrue(e_src, "TeV");

        // Perform numerical integration by summing
        const int nE = 10000;
        double emin  = 1000.0;            //!< 1 GeV
        double emax  = 3.0 * etrue.MeV(); //!< Migra_max = 3.0
        double dE    = (emax-emin)/nE;
        double sum   = 0.0;
        double E_obs = emin;
        for (int i = 0; i < nE; ++i, E_obs += dE) {
            GEnergy ereco(E_obs,"MeV");
            double dp  = edisp(ereco, etrue, dir);
            sum       += dp * dE;
        }
        test_value(sum, 1.0, 0.001, "Response cube energy dispersion integration "
                                    " for "+etrue.print());
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test CTA background
 ***************************************************************************/
void TestGCTAModel::test_model_bgd(void)
{
    // Test void constuctor
    test_try("Test void constuctor");
    try {
        GCTAModelBackground model;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test XML constuctor
    test_try("Test XML constuctor");
    try {
        GXml xml(cta_bgd_rad_gauss_xml);
        const GXmlElement& lib = *xml.element("source_library", 0);
        const GXmlElement& src = *lib.element("source", 0);
        GCTAModelBackground model(src);
        test_value(model["Prefactor"].value(), 61.8e-6);
        test_value(model["Index"].value(), -1.85);
        test_value(model["PivotEnergy"].value(), 1.0e6);
        test_assert(model.is_constant(), "Model is expected to be constant.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test spectral constuctor
    /*
    test_try("Test spectral constuctor");
    try {
        GModelSpectralPlaw plaw(1.0, 0.0, GEnergy(1.0, "TeV"));
        GCTAModelBackground model(plaw);
        test_value(model["Prefactor"].value(), 1.0);
        test_value(model["Index"].value(), 0.0);
        test_value(model["PivotEnergy"].value(), 1.0e6);
        test_assert(model.is_constant(), "Model is expected to be constant.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    */

    // Test XML loading of cube background
    GModels models(cta_bgd_rad_gauss_xml);
    GModel* model = models["My model"];
    test_value((*model)["Prefactor"].value(), 61.8e-6);
    test_value((*model)["Index"].value(), -1.85);
    test_value((*model)["PivotEnergy"].value(), 1.0e6);
    test_assert(model->is_constant(), "Model is expected to be constant.");

    // Test XML saving and reloading of cube background
    models.save("test.xml");
    models.load("test.xml");
    model = models["My model"];
    test_value((*model)["Prefactor"].value(), 61.8e-6);
    test_value((*model)["Index"].value(), -1.85);
    test_value((*model)["PivotEnergy"].value(), 1.0e6);
    test_assert(model->is_constant(), "Model is expected to be constant.");

    // Return
    return;

}


/***********************************************************************//**
 * @brief Test CTA cube background
 ***************************************************************************/
void TestGCTAModel::test_model_cube_bgd(void)
{
    // Test void constuctor
    test_try("Test void constuctor");
    try {
        GCTAModelCubeBackground model;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test XML constuctor
    test_try("Test XML constuctor");
    try {
        GXml xml(cta_cube_bgd_xml);
        const GXmlElement& lib = *xml.element("source_library", 0);
        const GXmlElement& src = *lib.element("source", 0);
        GCTAModelCubeBackground model(src);
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
        GCTAModelCubeBackground model(plaw);
        test_value(model["Prefactor"].value(), 1.0);
        test_value(model["Index"].value(), 0.0);
        test_value(model["PivotEnergy"].value(), 1.0e6);
        test_assert(model.is_constant(), "Model is expected to be constant.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test XML loading of cube background
    GModels models(cta_cube_bgd_xml);
    GModel* model = models["CTABackgroundModel"];
    test_value((*model)["Prefactor"].value(), 1.0);
    test_value((*model)["Index"].value(), 0.0);
    test_value((*model)["PivotEnergy"].value(), 1.0e6);
    test_assert(model->is_constant(), "Model is expected to be constant.");

    // Test XML saving and reloading of cube background
    models.save("test.xml");
    models.load("test.xml");
    model = models["CTABackgroundModel"];
    test_value((*model)["Prefactor"].value(), 1.0);
    test_value((*model)["Index"].value(), 0.0);
    test_value((*model)["PivotEnergy"].value(), 1.0e6);
    test_assert(model->is_constant(), "Model is expected to be constant.");

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
 * @brief Test CTA IRF background model
 ***************************************************************************/
void TestGCTAModel::test_model_aeff_bgd(void)
{
    // Test void constuctor
    test_try("Test void constuctor");
    try {
        GCTAModelAeffBackground model;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test XML constuctor
    test_try("Test XML constuctor");
    try {
        GXml xml(cta_aeff_bgd_xml);
        const GXmlElement& lib = *xml.element("source_library", 0);
        const GXmlElement& src = *lib.element("source", 0);
        GCTAModelIrfBackground model(src);
        test_value(model["Prefactor"].value(), 1.0e-14);
        test_value(model["Index"].value(), -2.4);
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
        GModelSpectralPlaw plaw(1.0e-14, -2.4, GEnergy(1.0, "TeV"));
        GCTAModelAeffBackground model(plaw);
        test_value(model["Prefactor"].value(), 1.0e-14);
        test_value(model["Index"].value(), -2.4);
        test_value(model["PivotEnergy"].value(), 1.0e6);
        test_assert(model.is_constant(), "Model is expected to be constant.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test XML loading of instrumental background
    GModels models(cta_aeff_bgd_xml);
    GModel* model = models["My model"];
    test_value((*model)["Prefactor"].value(), 1e-14);
    test_value((*model)["Index"].value(), -2.4);
    test_value((*model)["PivotEnergy"].value(), 1.0e6);
    test_assert(model->is_constant(), "Model is expected to be constant.");

    // Test XML saving and reloading of instrumental background
    models.save("test.xml");
    models.load("test.xml");
    model = models["My model"];
    test_value((*model)["Prefactor"].value(), 1.0e-14);
    test_value((*model)["Index"].value(), -2.4);
    test_value((*model)["PivotEnergy"].value(), 1.0e6);
    test_assert(model->is_constant(), "Model is expected to be constant.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Checks handling of CTA event bin
 ***************************************************************************/
void TestGCTAObservation::test_event_bin(void)
{
    // Test event bin void constructor
    GCTAEventBin bin;
    test_value(bin.classname(), "GCTAEventBin",
               "Check classname() for empty bin");
    test_value(bin.size(), 0.0, 1.0e-10, "Check size() for empty bin");
    test_assert(!bin.dir().has_dir(), "Check dir().has_dir() for empty bin");
    test_assert(!bin.dir().has_detx(), "Check dir().has_detx() for empty bin");
    test_assert(!bin.dir().has_dety(), "Check dir().has_dety() for empty bin");
    test_value(bin.energy().MeV(), 0.0, 1.0e-10,
               "Check energy() for empty bin");
    test_value(bin.time().secs(), 0.0, 1.0e-10, "Check time() for empty bin");
    test_value(bin.counts(), 0.0, 1.0e-10, "Check counts() for empty bin");
    test_value(bin.error(),  0.0, 1.0e-10, "Check error() for empty bin");
    test_value(bin.ipix(), -1, "Check ipix() for empty bin");
    test_value(bin.ieng(), -1, "Check ieng() for empty bin");
    test_value(bin.solidangle(), 0.0, 1.0e-10,
               "Check solidangle() for empty bin");
    test_value(bin.ewidth().MeV(), 0.0, 1.0e-10,
               "Check ewidth() for empty bin");
    test_value(bin.ontime(), 0.0, 1.0e-10, "Check ontime() for empty bin");
    test_value(bin.weight(), 0.0, 1.0e-10, "Check weight() for empty bin");
    test_value(bin.print(), "0", "Check print() for empty bin");

    // Set bin attributes
    GSkyDir  skydir;
    skydir.radec_deg(37.2, -67.3);
    GCTAInstDir instdir(skydir);
    instdir.detx(-0.78);
    instdir.dety(+0.35);
    bin.dir(instdir);
    bin.energy(GEnergy(2.0, "TeV"));
    bin.time(GTime(87.3));
    bin.counts(4.0);
    bin.ipix(5);
    bin.ieng(7);
    bin.solidangle(3.14);
    bin.ewidth(GEnergy(7.0, "GeV"));
    bin.ontime(101.0);
    bin.weight(0.71);

    // Set size reference
    double ref_size = 3.14 * 7000.0 * 101.0 * 0.71;

    // Test bin attributes
    test_value(bin.size(), ref_size, 1.0e-10, "Check size() for filled bin");
    test_assert(bin.dir().has_dir(), "Check dir().has_dir() for filled bin");
    test_assert(bin.dir().has_detx(), "Check dir().has_detx() for filled bin");
    test_assert(bin.dir().has_dety(), "Check dir().has_dety() for filled bin");
    test_value(bin.dir().dir().ra_deg(), 37.2, 1.0e-10,
               "Check dir().dir().ra_deg() for filled bin");
    test_value(bin.dir().dir().dec_deg(), -67.3, 1.0e-10,
               "Check dir().dir().dec_deg() for filled bin");
    test_value(bin.dir().detx(), -0.78, 1.0e-10,
               "Check dir().detx() for filled bin");
    test_value(bin.dir().dety(), +0.35, 1.0e-10,
               "Check dir().detx() for filled bin");
    test_value(bin.energy().TeV(), 2.0, 1.0e-10,
               "Check energy() for filled bin");
    test_value(bin.time().secs(), 87.3, 1.0e-10, "Check time() for filled bin");
    test_value(bin.counts(), 4.0, 1.0e-10, "Check counts() for filled bin");
    test_value(bin.error(),  2.0, 1.0e-10, "Check error() for filled bin");
    test_value(bin.ipix(), 5, "Check ipix() for filled bin");
    test_value(bin.ieng(), 7, "Check ieng() for filled bin");
    test_value(bin.solidangle(), 3.14, 1.0e-10,
               "Check solidangle() for filled bin");
    test_value(bin.ewidth().GeV(), 7.0, 1.0e-10,
               "Check ewidth() for filled bin");
    test_value(bin.ontime(), 101.0, 1.0e-10, "Check ontime() for filled bin");
    test_value(bin.weight(), 0.71, 1.0e-10, "Check weight() for filled bin");
    test_value(bin.print(), "4", "Check print() for filled bin");

    // Test copy constructor
    GCTAEventBin bin2(bin);
    test_value(bin2.size(), ref_size, 1.0e-10, "Check size() for copied bin");
    test_value(bin2.dir().dir().ra_deg(), 37.2, 1.0e-10,
               "Check dir().dir().ra_deg() for filled bin");
    test_value(bin2.dir().dir().dec_deg(), -67.3, 1.0e-10,
               "Check dir().dir().dec_deg() for filled bin");
    test_value(bin2.dir().detx(), -0.78, 1.0e-10,
               "Check dir().detx() for filled bin");
    test_value(bin2.dir().dety(), +0.35, 1.0e-10,
               "Check dir().detx() for filled bin");
    test_value(bin2.energy().TeV(), 2.0, 1.0e-10,
               "Check energy() for copied bin");
    test_value(bin2.time().secs(), 87.3, 1.0e-10, "Check time() for copied bin");
    test_value(bin2.counts(), 4.0, 1.0e-10, "Check counts() for copied bin");
    test_value(bin2.error(),  2.0, 1.0e-10, "Check error() for copied bin");
    test_value(bin2.ipix(), 5, "Check ipix() for copied bin");
    test_value(bin2.ieng(), 7, "Check ieng() for copied bin");
    test_value(bin2.solidangle(), 3.14, 1.0e-10,
               "Check solidangle() for copied bin");
    test_value(bin2.ewidth().GeV(), 7.0, 1.0e-10,
               "Check ewidth() for copied bin");
    test_value(bin2.ontime(), 101.0, 1.0e-10, "Check ontime() for copied bin");
    test_value(bin2.weight(), 0.71, 1.0e-10, "Check weight() for copied bin");
    test_value(bin2.print(), "4", "Check print() for copied bin");

    // Assignment operator
    GCTAEventBin bin3 = bin;
    test_value(bin3.size(), ref_size, 1.0e-10, "Check size() for assigned bin");
    test_value(bin3.dir().dir().ra_deg(), 37.2, 1.0e-10,
               "Check dir().dir().ra_deg() for filled bin");
    test_value(bin3.dir().dir().dec_deg(), -67.3, 1.0e-10,
               "Check dir().dir().dec_deg() for filled bin");
    test_value(bin3.dir().detx(), -0.78, 1.0e-10,
               "Check dir().detx() for filled bin");
    test_value(bin3.dir().dety(), +0.35, 1.0e-10,
               "Check dir().detx() for filled bin");
    test_value(bin3.energy().TeV(), 2.0, 1.0e-10,
               "Check energy() for assigned bin");
    test_value(bin3.time().secs(), 87.3, 1.0e-10,
               "Check time() for assigned bin");
    test_value(bin3.counts(), 4.0, 1.0e-10, "Check counts() for assigned bin");
    test_value(bin3.error(),  2.0, 1.0e-10, "Check error() for assigned bin");
    test_value(bin3.ipix(), 5, "Check ipix() for assigned bin");
    test_value(bin3.ieng(), 7, "Check ieng() for assigned bin");
    test_value(bin3.solidangle(), 3.14, 1.0e-10,
               "Check solidangle() for assigned bin");
    test_value(bin3.ewidth().GeV(), 7.0, 1.0e-10,
               "Check ewidth() for assigned bin");
    test_value(bin3.ontime(), 101.0, 1.0e-10, "Check ontime() for assigned bin");
    test_value(bin3.weight(), 0.71, 1.0e-10, "Check weight() for assigned bin");
    test_value(bin3.print(), "4", "Check print() for assigned bin");

    // clone method
    GCTAEventBin* bin4 = bin.clone();
    test_value(bin4->size(), ref_size, 1.0e-10, "Check size() for cloned bin");
    test_value(bin4->dir().dir().ra_deg(), 37.2, 1.0e-10,
               "Check dir().dir().ra_deg() for filled bin");
    test_value(bin4->dir().dir().dec_deg(), -67.3, 1.0e-10,
               "Check dir().dir().dec_deg() for filled bin");
    test_value(bin4->dir().detx(), -0.78, 1.0e-10,
               "Check dir().detx() for filled bin");
    test_value(bin4->dir().dety(), +0.35, 1.0e-10,
               "Check dir().detx() for filled bin");
    test_value(bin4->energy().TeV(), 2.0, 1.0e-10,
               "Check energy() for cloned bin");
    test_value(bin4->time().secs(), 87.3, 1.0e-10,
               "Check time() for cloned bin");
    test_value(bin4->counts(), 4.0, 1.0e-10, "Check counts() for cloned bin");
    test_value(bin4->error(),  2.0, 1.0e-10, "Check error() for cloned bin");
    test_value(bin4->ipix(), 5, "Check ipix() for cloned bin");
    test_value(bin4->ieng(), 7, "Check ieng() for cloned bin");
    test_value(bin4->solidangle(), 3.14, 1.0e-10,
               "Check solidangle() for cloned bin");
    test_value(bin4->ewidth().GeV(), 7.0, 1.0e-10,
               "Check ewidth() for cloned bin");
    test_value(bin4->ontime(), 101.0, 1.0e-10, "Check ontime() for cloned bin");
    test_value(bin4->weight(), 0.71, 1.0e-10, "Check weight() for cloned bin");
    test_value(bin4->print(), "4", "Check print() for cloned bin");

    // clear method
    bin.clear();
    test_value(bin.classname(), "GCTAEventBin",
               "Check classname() for cleared bin");
    test_value(bin.size(), 0.0, 1.0e-10, "Check size() for cleared bin");
    test_assert(!bin.dir().has_dir(), "Check dir().has_dir() for cleared bin");
    test_assert(!bin.dir().has_detx(), "Check dir().has_detx() for cleared bin");
    test_assert(!bin.dir().has_dety(), "Check dir().has_dety() for cleared bin");
    test_value(bin.energy().MeV(), 0.0, 1.0e-10,
               "Check energy() for cleared bin");
    test_value(bin.time().secs(), 0.0, 1.0e-10, "Check time() for cleared bin");
    test_value(bin.counts(), 0.0, 1.0e-10, "Check counts() for cleared bin");
    test_value(bin.error(),  0.0, 1.0e-10, "Check error() for cleared bin");
    test_value(bin.ipix(), -1, "Check ipix() for cleared bin");
    test_value(bin.ieng(), -1, "Check ieng() for cleared bin");
    test_value(bin.solidangle(), 0.0, 1.0e-10,
               "Check solidangle() for cleared bin");
    test_value(bin.ewidth().MeV(), 0.0, 1.0e-10,
               "Check ewidth() for cleared bin");
    test_value(bin.ontime(), 0.0, 1.0e-10, "Check ontime() for cleared bin");
    test_value(bin.weight(), 0.0, 1.0e-10, "Check weight() for cleared bin");
    test_value(bin.print(), "0", "Check print() for cleared bin");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test event cube
 ***************************************************************************/
void TestGCTAObservation::test_event_cube(void)
{
    // Construct empty event cube
    GCTAEventCube cube1;

    // Test empty event cube
    test_value(cube1.size(), 0, "Check that empty event cube has size()=0");
    test_value(cube1.dim(), 0, "Check that empty event cube has dim()=0");
    test_value(cube1.number(), 0, "Check that empty event cube has number()=0");
    test_value(cube1.nx(), 0, "Check that empty event cube has nx()=0");
    test_value(cube1.ny(), 0, "Check that empty event cube has ny()=0");
    test_value(cube1.npix(), 0, "Check that empty event cube has npix()=0");
    test_value(cube1.ebins(), 0, "Check that empty event cube has ebins()=0");

    // Construct 3x3x5 event cube
    GSkyMap       map("CAR", "CEL", 0.0, 0.0, 0.5, 0.5, 3, 4, 5);
    GEbounds      ebounds(5, GEnergy(0.1, "TeV"), GEnergy(10.0, "TeV"));
    GGti          gti(GTime(0.0, "sec"), GTime(1800.0, "sec"));
    GCTAEventCube cube2(map, ebounds, gti);

    // Test 3x3x5 event cube
    test_value(cube2.size(), 60, "Check that event cube has size()=60");
    test_value(cube2.dim(), 3, "Check that event cube has dim()=3");
    test_value(cube2.number(), 0, "Check that event cube has number()=0");
    test_value(cube2.nx(), 3, "Check that event cube has nx()=3");
    test_value(cube2.ny(), 4, "Check that event cube has ny()=4");
    test_value(cube2.npix(), 12, "Check that event cube has npix()=12");
    test_value(cube2.ebins(), 5, "Check that event cube has ebins()=5");

    // Fill counts cube and test result
    double ref(0.0);
    double value(0.0);
    for (int ix = 0; ix < cube2.nx(); ++ix) {
        for (int iy = 0; iy < cube2.ny(); ++iy) {
            for (int iebin = 0; iebin < cube2.ebins(); ++iebin) {
                value += 0.376;
                ref   += value;
                map(GSkyPixel(ix,iy),iebin) = value;
            }
        }
    }
    cube2.counts(map);
    test_value(cube2.number(), int(ref+0.5),
               "Check event cube number() for filled cube");

    // Checks weights
    double total(0.0);
    for (int ix = 0; ix < cube2.nx(); ++ix) {
        for (int iy = 0; iy < cube2.ny(); ++iy) {
            for (int iebin = 0; iebin < cube2.ebins(); ++iebin) {
                total += cube2.weights()(GSkyPixel(ix,iy),iebin);
            }
        }
    }
    test_value(total, 60.0, 1.0e-6,
               "Check event cube weights for cube with unity weights");

    // Save event cube
    cube2.save("test_cta_event_cube.fits", true);

    // Construct counts cube from FITS file
    GCTAEventCube cube3("test_cta_event_cube.fits");

    // Test loaded event cube
    test_value(cube3.size(), 60, "Check that loaded event cube has size()=60");
    test_value(cube3.dim(), 3, "Check that loaded event cube has dim()=3");
    test_value(cube3.nx(), 3, "Check that loaded event cube has nx()=3");
    test_value(cube3.ny(), 4, "Check that loaded event cube has ny()=4");
    test_value(cube3.npix(), 12, "Check that loaded event cube has npix()=12");
    test_value(cube3.ebins(), 5, "Check that loaded event cube has ebins()=5");
    test_value(cube3.number(), int(ref+0.5),
               "Check event cube number() for loaded cube");

    // Checks weights for loaded event cube
    total = 0.0;
    for (int ix = 0; ix < cube2.nx(); ++ix) {
        for (int iy = 0; iy < cube2.ny(); ++iy) {
            for (int iebin = 0; iebin < cube2.ebins(); ++iebin) {
                total += cube3.weights()(GSkyPixel(ix,iy),iebin);
            }
        }
    }
    test_value(total, 60.0, 1.0e-6,
               "Check event cube weights for loaded cube with unity weights");

    // Set weights
    cube3.weights(map);

    // Checks weights
    total = 0.0;
    for (int ix = 0; ix < cube2.nx(); ++ix) {
        for (int iy = 0; iy < cube2.ny(); ++iy) {
            for (int iebin = 0; iebin < cube2.ebins(); ++iebin) {
                total += cube3.weights()(GSkyPixel(ix,iy),iebin);
            }
        }
    }
    test_value(total, ref, 1.0e-6,
               "Check event cube weights after setting them");

    // Save event cube
    cube3.save("test_cta_event_cube.fits", true);

    // Construct counts cube from FITS file
    GCTAEventCube cube4("test_cta_event_cube.fits");

    // Checks weights
    total = 0.0;
    for (int ix = 0; ix < cube2.nx(); ++ix) {
        for (int iy = 0; iy < cube2.ny(); ++iy) {
            for (int iebin = 0; iebin < cube2.ebins(); ++iebin) {
                total += cube3.weights()(GSkyPixel(ix,iy),iebin);
            }
        }
    }
    test_value(total, ref, 1.0e-6,
               "Check event cube weights after loading them again");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test unbinned observation handling
 ***************************************************************************/
void TestGCTAObservation::test_unbinned_obs(void)
{
    // Declare observations
    GObservations   obs;
    GCTAObservation run;

    // Load unbinned CTA observation
    test_try("Load unbinned CTA observation");
    try {
        run.load(cta_events);
        run.response(cta_irf, GCaldb(cta_caldb));
        test_value(run.roi().centre().dir().ra_deg(), 83.63);
        test_value(run.roi().centre().dir().dec_deg(), 22.01);
        test_value(run.roi().radius(), 5.0);
        test_value(run.ebounds().emin().TeV(), 0.1);
        test_value(run.ebounds().emax().TeV(), 100.0);
        test_value(run.gti().tstart().convert(run.gti().reference()), 0.0);
        test_value(run.gti().tstop().convert(run.gti().reference()), 1800.0);
        test_value(run.ontime(), 1800.0);
        test_value(run.livetime(), 1710.0);
        test_value(run.deadc(), 0.95);
        test_value(run.ra_obj(), 0.0);
        test_value(run.dec_obj(), 0.0);
        test_value(run.id(), "0");
        test_value(run.pointing().dir().ra_deg(), 83.63);
        test_value(run.pointing().dir().dec_deg(), 22.01);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Load unbinned CTA observation with non-default extension names
    test_try("Load unbinned CTA observation with non-default extension names");
    try {
        run.load(cta_events_gti);
        run.response(cta_irf, GCaldb(cta_caldb));
        test_value(run.roi().centre().dir().ra_deg(), 83.65);
        test_value(run.roi().centre().dir().dec_deg(), 23.01);
        test_value(run.roi().radius(), 4.0);
        test_value(run.ebounds().emin().TeV(), 0.2);
        test_value(run.ebounds().emax().TeV(), 120.0);
        test_value(run.gti().tstart().convert(run.gti().reference()), 1.0);
        test_value(run.gti().tstop().convert(run.gti().reference()), 2000.0);
        test_value(run.ontime(), 1800.0);
        test_value(run.livetime(), 1710.0);
        test_value(run.deadc(), 0.95);
        test_value(run.ra_obj(), 0.0);
        test_value(run.dec_obj(), 0.0);
        test_value(run.id(), "0");
        test_value(run.pointing().dir().ra_deg(), 83.63);
        test_value(run.pointing().dir().dec_deg(), 22.01);
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
    test_value(num, 4134, 1.0e-20, "Test event iterator");

    // Test XML loading and saving
    test_try("Test XML loading and saving");
    try {
        obs = GObservations(cta_unbin_xml);
        obs.save("test_cta_obs_unbinned.xml");
        GCTAObservation* run = dynamic_cast<GCTAObservation*>(obs[0]);
        test_value(run->roi().centre().dir().ra_deg(), 83.63);
        test_value(run->roi().centre().dir().dec_deg(), 22.01);
        test_value(run->roi().radius(), 5.0);
        test_value(run->ebounds().emin().TeV(), 0.1);
        test_value(run->ebounds().emax().TeV(), 100.0);
        test_value(run->gti().tstart().convert(run->gti().reference()), 0.0);
        test_value(run->gti().tstop().convert(run->gti().reference()), 1800.0);
        test_value(run->ontime(), 1800.0);
        test_value(run->livetime(), 1710.0);
        test_value(run->deadc(), 0.95);
        test_value(run->ra_obj(), 0.0);
        test_value(run->dec_obj(), 0.0);
        test_value(run->id(), "00001");
        test_value(run->pointing().dir().ra_deg(), 83.63);
        test_value(run->pointing().dir().dec_deg(), 22.01);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Remove test files
    std::remove("test_cta_events1.fits");
    std::remove("test_cta_events2.fits");

    // Test loading of the event file and saving into a different extension
    run.load(cta_events);
    run.save("test_cta_events1.fits[EVENTS2;GTI2]", true);
    GFits fits("test_cta_events1.fits");
    test_value(fits.size(), 3, "FITS file should contain 3 HDUs");
    test_assert(fits.contains("EVENTS2"), "FITS should contain \"EVENTS2\" HDU");
    test_assert(fits.contains("GTI2"), "FITS should contain \"GTI2\" HDU");
    fits.close();
    run.save("test_cta_events1.fits[EVENTS3;GTI3]", true);
    fits.open("test_cta_events1.fits");
    test_value(fits.size(), 5, "FITS file should contain 5 HDUs");
    test_assert(fits.contains("EVENTS2"), "FITS should contain \"EVENTS2\" HDU");
    test_assert(fits.contains("GTI2"), "FITS should contain \"GTI2\" HDU");
    test_assert(fits.contains("EVENTS3"), "FITS should contain \"EVENTS3\" HDU");
    test_assert(fits.contains("GTI3"), "FITS should contain \"GTI3\" HDU");
    fits.close();

    // Test writing of multiple event files into FITS file
    run.load(cta_events);
    GFits fits2;
    run.write(fits2, "EVENTS1", "GTI1");
    run.write(fits2, "EVENTS2", "GTI2");
    run.write(fits2, "EVENTS3", "GTI3");
    fits2.saveto("test_cta_events2.fits", true);
    fits2.close();
    fits.open("test_cta_events2.fits");
    test_value(fits.size(), 7, "FITS file should contain 7 HDUs");
    test_assert(fits.contains("EVENTS1"), "FITS should contain \"EVENTS1\" HDU");
    test_assert(fits.contains("GTI1"), "FITS should contain \"GTI1\" HDU");
    test_assert(fits.contains("EVENTS2"), "FITS should contain \"EVENTS2\" HDU");
    test_assert(fits.contains("GTI2"), "FITS should contain \"GTI2\" HDU");
    test_assert(fits.contains("EVENTS3"), "FITS should contain \"EVENTS3\" HDU");
    test_assert(fits.contains("GTI3"), "FITS should contain \"GTI3\" HDU");
    fits.close();

    // Return
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

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test stacked observation handling
 ***************************************************************************/
void TestGCTAObservation::test_stacked_obs(void)
{
    // Construct stacked observation without energy dispersion
    GCTAObservation cta1(cta_stacked_cntcube,
                         cta_stacked_expcube,
                         cta_stacked_psfcube,
                         cta_stacked_bkgcube);

    // Test for presence of response
    const GCTAResponseCube* rsp =
          dynamic_cast<const GCTAResponseCube*>(cta1.response());
    test_assert((rsp != NULL), "Observation contains cube response");
    if (rsp != NULL) {
        rsp->apply_edisp(true);  // Try to use energy dispersion
        test_assert((!rsp->use_edisp()), "Response has no energy dispersion");
    }

    // Construct stacked observation with energy dispersion
    GCTAObservation cta2(cta_stacked_cntcube,
                         cta_stacked_expcube,
                         cta_stacked_psfcube,
                         cta_stacked_edispcube,
                         cta_stacked_bkgcube);

    // Test for presence of response
    rsp = dynamic_cast<const GCTAResponseCube*>(cta2.response());
    test_assert((rsp != NULL), "Observation contains cube response");
    if (rsp != NULL) {
        rsp->apply_edisp(true);  // Try to use energy dispersion
        test_assert(rsp->use_edisp(), "Response has energy dispersion");
    }

    // Construct stacked observation from XML file without energy dispersion
    GObservations obs(cta_stacked_xml);
    test_value(obs.size(), 1, "One observation in container");

    // Test for presence of response
    rsp = dynamic_cast<const GCTAResponseCube*>(obs[0]->response());
    test_assert((rsp != NULL), "Observation contains cube response");
    if (rsp != NULL) {
        rsp->apply_edisp(true);  // Try to use energy dispersion
        test_assert((!rsp->use_edisp()), "Response has no energy dispersion");
    }

    // Save observation container into XML file
    obs.save("test_cta_obs_cube.xml");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test On/Off observation handling
 ***************************************************************************/
void TestGCTAObservation::test_onoff_obs(void)
{
    // Load On/Off observation into container
    GObservations obs(cta_onoff_obs);

    // Save observation container into XML file
    obs.save("test_cta_onoff_obs.xml");

    // Load CTA observations
    obs.load(cta_unbin_xml);

    // Setup sky regions
    GSkyRegions regions(cta_onoff_onreg);

    // Set and check regions
    for (int i = 0; i < obs.size(); ++i) {
        static_cast<GCTAObservation*>(obs[i])->off_regions(regions);
        test_value(static_cast<const GCTAObservation*>(obs[i])->off_regions().size(),
                   1, "Check number of sky regions");
    }

    // Save observation container into XML file
    obs.save("test_cta_onoff_obs_regions.xml");

    // Re-load observation container and check regions
    obs.load("test_cta_onoff_obs_regions.xml");
    
    // Check regions
    for (int i = 0; i < obs.size(); ++i) {
        test_value(static_cast<const GCTAObservation*>(obs[i])->off_regions().size(),
                   1, "Check number of sky regions");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test unbinned optimizer
 ***************************************************************************/
void TestGCTAOptimize::test_unbinned_optimizer(void)
{
    // Set reference result
    double fit_results[] = {83.6331, 0,
                            22.0145, 0,
                            6.03743e-16, 2.01634e-17,
                            -2.49595, 0.0249533,
                            300000, 0,
                            1, 0,
                            2.94185, 0.0379307,
                            6.43398e-05, 1.80715e-06,
                            -1.82866, 0.0155452,
                            1000000, 0,
                            1, 0};

    // Declare observations
    GObservations   obs;
    GCTAObservation run;

    // Load unbinned CTA observation and set response
    run.load(cta_events);
    run.response(cta_irf, GCaldb(cta_caldb));
    obs.append(run);

    // Load models from XML file
    obs.models(cta_model_xml);

    // Perform LM optimization
    GOptimizerLM opt;
    opt.max_iter(100);
    obs.optimize(opt);
    obs.errors(opt);

    // Verify fit results
    check_results(obs, fit_results);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test binned optimizer
 ***************************************************************************/
void TestGCTAOptimize::test_binned_optimizer(void)
{
    // Set reference result
    double fit_results[] = {83.6331, 0,
                            22.0145, 0,
                            5.99255e-16, 2.00626e-17,
                            -2.49175, 0.0250639,
                            300000, 0,
                            1, 0,
                            2.95662, 0.0704027,
                            6.42985e-05, 1.96478e-06,
                            -1.82084, 0.0163749,
                            1000000, 0,
                            1, 0};

    // Declare observations
    GObservations   obs;
    GCTAObservation run;

    // Load binned CTA observation and set response
    run.load(cta_cntmap);
    run.response(cta_irf, GCaldb(cta_caldb));
    obs.append(run);

    // Load models from XML file
    obs.models(cta_model_xml);

    // Perform LM optimization
    GOptimizerLM opt;
    opt.max_iter(100);
    obs.optimize(opt);
    obs.errors(opt);

    // Verify fit results
    check_results(obs, fit_results);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test stacked optimizer
 *
 * The stacked response cubes have been computed using the Prod2::South_0.5h
 * response.
 ***************************************************************************/
void TestGCTAOptimize::test_stacked_optimizer(void)
{
    // Set reference result
    double fit_results[] = {83.6331, 0,
                            22.0145, 0,
                            5.98707105445563e-16, 1.09127928921077e-17,
                            -2.50267020665246, 0.0145532902844396,
                            300000, 0,
                            1, 0,
                            1.61294446311705, 0.033167560400167,
                            -0.218565245795573, 0.0146770868941243,
                            1.0e6, 0,
                            1, 0};

    // Load stacked CTA observation
    GObservations obs(cta_stacked_xml);

    // Load models from XML file
    obs.models(cta_stacked_model);

    // Perform LM optimization
    GOptimizerLM opt;
    opt.max_iter(100);
    obs.optimize(opt);
    obs.errors(opt);

    // Verify fit results
    check_results(obs, fit_results);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test On/Off optimizer
 ***************************************************************************/
void TestGCTAOptimize::test_onoff_optimizer_cstat(void)
{
    // Reference result
    double fit_results[] = {83.6331, 0,
                            22.0145, 0,
                            5.745479e-16, 3.061254e-17,
                            -2.488965, 0.021807,
                            300000, 0,
                            1, 0,
                            1.002526, 0.061438,
                            -0.044025, 0.055257,
                            1.0e6, 0,
                            1, 0};

    // Load On/Off CTA observation
    GObservations obs(cta_onoff_obs);

    // Set fit statistic
    for (int i = 0; i < obs.size(); ++i){
        obs[i]->statistic("CSTAT");
    }

    // Load models from XML file
    obs.models(cta_onoff_model);

    // Perform LM optimization
    GOptimizerLM opt;
    opt.max_iter(100);
    obs.optimize(opt);
    obs.errors(opt);

    // Verify fit results
    check_results(obs, fit_results);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test On/Off optimizer
 ***************************************************************************/
void TestGCTAOptimize::test_onoff_optimizer_wstat(void)
{
    // Reference result
    double fit_results[] = {83.6331, 0,
                            22.0145, 0,
                            5.746472e-16, 3.063591e-17,
                            -2.488890, 0.021824,
                            300000, 0,
                            1, 0,
                            1.0, 0.,
                            0., 0.,
                            1.0e6, 0,
                            1, 0};

    // Load On/Off CTA observation
    GObservations obs(cta_onoff_obs);

    // Set fit statistic
    for (int i = 0; i < obs.size(); ++i){
        obs[i]->statistic("WSTAT");
    }

    // Load models from XML file
    obs.models(cta_onoff_model);

    // Perform LM optimization
    GOptimizerLM opt;
    opt.max_iter(100);
    obs.optimize(opt);
    obs.errors(opt);

    // Verify fit results
    check_results(obs, fit_results);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Check fit results
 *
 * @param[in] obs Observation container
 * @param[in] results Pointer to array of expected fit results
 ***************************************************************************/
void TestGCTAOptimize::check_results(const GObservations& obs,
                                     const double*        results)
{
    // Verify fit results
    for (int i = 0, j = 0; i < obs.models().size(); ++i) {
        const GModel* model = obs.models()[i];
        for (int k = 0; k < model->size(); ++k) {
            GModelPar par  = (*model)[k];
            std::string msg = "Verify optimization result for " + par.print();
            test_value(par.value(), results[j],
                       std::abs(1.0e-2*results[j]), msg);
            j++;
            test_value(par.error(), results[j],
                       std::abs(1.0e-2*results[j]), msg);
            j++;
        }
    }

    // Return
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

    // Set CALDB environment variable
    setenv("CALDB", cta_caldb.c_str(), 1);

    // Initially assume that we pass all tests
    bool success = true;

    // Create test suites and append them to the container
    TestGCTA                test;
    TestGCTAResponse        rsp;
    TestGCTAModel           model;
    TestGCTAOptimize        opt;
    TestGCTAObservation     obs;
    testsuites.append(rsp);
    if (has_data) {
        testsuites.append(test);
    	testsuites.append(model);
        testsuites.append(opt);
        testsuites.append(obs);
    }

    // Run the testsuites
    success = testsuites.run();

    // Save test report
    testsuites.save("reports/GCTA.xml");

    // Return success status
    return (success ? 0 : 1);
}
