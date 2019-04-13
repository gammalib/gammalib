/***************************************************************************
 *                    test_CTA.hpp  -  Test CTA classes                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2019 by Juergen Knoedlseder                         *
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
 * @file test_CTA.hpp
 * @brief Definition of CTA classes unit tests
 * @author Juergen Knoedlseder
 */

#ifndef TEST_CTA_HPP
#define TEST_CTA_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"


/***********************************************************************//**
 * @class TestGCTA
 *
 * @brief Test suite for miscellaneous CTA classes
 *
 * This class defines a unit test suite for miscellaneous CTA classes
 ***************************************************************************/
class TestGCTA : public GTestSuite {
public:

    // Constructors and destructors
    TestGCTA(void) : GTestSuite() {}
    virtual ~TestGCTA(void) {}

    // Methods
    virtual void        set(void);
    virtual TestGCTA*   clone(void) const { return new TestGCTA(*this); }
    virtual std::string classname(void) const { return "TestGCTA"; }
    void                test_instdir(void);
    void                test_pointing(void);
};


/***********************************************************************//**
 * @class TestGCTAResponse
 *
 * @brief Test suite for CTA response class testing
 *
 * This class defines a unit test suite for the CTA response class.
 ***************************************************************************/
class TestGCTAResponse : public GTestSuite {
public:
    // Constructors and destructors
    TestGCTAResponse(void) : GTestSuite() {}
    virtual ~TestGCTAResponse(void) {}

    // Methods
    virtual void              set(void);
    virtual TestGCTAResponse* clone(void) const { return new TestGCTAResponse(*this); }
    virtual std::string       classname(void) const { return "TestGCTAResponse"; }
    void                      test_response(void);
    void                      test_response_aeff(void);
    void                      test_response_psf(void);
    void                      test_response_psf_king(void);
    void                      test_response_psf_table(void);
    void                      test_response_npsf(void);
    void                      test_response_irf_diffuse(void);
    void                      test_response_npred_diffuse(void);
    void                      test_response_edisp(void);
    void                      test_response_edisp_PerfTable(void);
    void                      test_response_edisp_RMF(void);
    void                      test_response_edisp_2D(void);
    void                      test_response_bgd_PerfTable(void);
    void                      test_response_bgd_3D(void);
    void                      test_response_expcube(void);
    void                      test_response_psfcube(void);
    void                      test_response_bkgcube(void);
    void                      test_response_edispcube(void);
    void                      test_response_cache(void);

    // Utility methods
    void test_response_edisp_integration(const GCTAResponseIrf& rsp,
                                         const double&          e_src_min = 0.1,
                                         const double&          e_src_max = 10.0);
    void test_edisp_integration(const GCTAEdisp&   edisp,
                                const double&      e_src_min = 0.1,
                                const double&      e_src_max = 10.0);
    void test_edispcube_integration(const GCTACubeEdisp& edisp,
                                    const double&        e_src_min = 0.1,
                                    const double&        e_src_max = 10.0);
};


/***********************************************************************//**
 * @class TestGCTAModel
 *
 * @brief Test suite for CTA model class testing
 *
 * This class defines a unit test suite for the CTA model class.
 ***************************************************************************/
class TestGCTAModel : public GTestSuite {
public:
    // Constructors and destructors
	TestGCTAModel(void) : GTestSuite() {}
    virtual ~TestGCTAModel(void) {}

    // Methods
    virtual void           set(void);
    virtual TestGCTAModel* clone(void) const { return new TestGCTAModel(*this); }
    virtual std::string    classname(void) const { return "TestGCTAModel"; }
    void                   test_model_bgd(void);
    void                   test_model_cube_bgd(void);
    void                   test_model_irf_bgd(void);
    void                   test_model_aeff_bgd(void);
    void                   test_spatial_gauss_spectrum(void);
};


/***********************************************************************//**
 * @class TestGCTAObservation
 *
 * @brief Test suite for CTA observation testing
 *
 * This class defines a unit test suite for testing of GCTAObservation.
 ***************************************************************************/
class TestGCTAObservation : public GTestSuite {
public:
    // Constructors and destructors
    TestGCTAObservation(void) : GTestSuite() {}
    virtual ~TestGCTAObservation(void) {}

    // Methods
    virtual void                 set(void);
    virtual TestGCTAObservation* clone(void) const { return new TestGCTAObservation(*this); }
    virtual std::string          classname(void) const { return "TestGCTAObservation"; }
    void                         test_event_bin(void);
    void                         test_event_cube(void);
    void                         test_unbinned_obs(void);
    void                         test_binned_obs(void);
    void                         test_stacked_obs(void);
    void                         test_onoff_obs(void);
};


/***********************************************************************//**
 * @class TestGCTAOptimize
 *
 * @brief Test suite for CTA optimizer testing
 *
 * This class defines a unit test suite for testing of CTA optimizers.
 ***************************************************************************/
class TestGCTAOptimize : public GTestSuite {
public:
    // Constructors and destructors
    TestGCTAOptimize(void) : GTestSuite() {}
    virtual ~TestGCTAOptimize(void) {}

    // Methods
    virtual void              set(void);
    virtual TestGCTAOptimize* clone(void) const { return new TestGCTAOptimize(*this); }
    virtual std::string       classname(void) const { return "TestGCTAOptimize"; }
    void                      test_unbinned_optimizer(void);
    void                      test_binned_optimizer(void);
    void                      test_stacked_optimizer(void);
    void                      test_onoff_optimizer_cstat(void);
    void                      test_onoff_optimizer_wstat(void);

protected:
    // Protected methods
    void check_results(const GObservations& obs, const double* results);
};

#endif /* TEST_CTA_HPP */
