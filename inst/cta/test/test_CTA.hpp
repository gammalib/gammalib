/***************************************************************************
 *                    test_CTA.hpp  -  Test CTA classes                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2015 by Juergen Knoedlseder                         *
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
    virtual TestGCTAResponse* clone(void) const;
    virtual std::string       classname(void) const { return "TestGCTAResponse"; }
    void                      test_response(void);
    void                      test_response_aeff(void);
    void                      test_response_psf(void);
    void                      test_response_psf_king(void);
    void                      test_response_npsf(void);
    void                      test_response_edisp(void);
    void                      test_response_edisp_PerfTable(void);
    void                      test_response_edisp_RMF(void);
    void                      test_response_edisp_2D(void);
    void                      test_response_irf_diffuse(void);
    void                      test_response_npred_diffuse(void);
    void                      test_response_expcube(void);
    void                      test_response_psfcube(void);
    void                      test_response_bkgcube(void);
    void                      test_response_edispcube(void);

    // Utility methods
    void test_response_edisp_integration(const GCTAResponseIrf& rsp,
                                         const double&          e_src_min = 0.1,
                                         const double&          e_src_max = 10.0);
    void test_edisp_integration(const GCTAEdisp& edisp,
                                const double&    e_src_min = 0.1,
                                const double&    e_src_max = 10.0);
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
    virtual TestGCTAModel* clone(void) const;
    virtual std::string    classname(void) const { return "TestGCTAModel"; }
    void                   test_model_cube_bgd(void);
    void                   test_model_irf_bgd(void);
    void                   test_model_aeff_bgd(void);
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
    virtual TestGCTAObservation* clone(void) const;
    virtual std::string          classname(void) const { return "TestGCTAObservation"; }
    void                         test_unbinned_obs(void);
    void                         test_binned_obs(void);
    void                         test_stacked_obs(void);
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
    virtual TestGCTAOptimize* clone(void) const;
    virtual std::string       classname(void) const { return "TestGCTAOptimize"; }
    void                      test_unbinned_optimizer(void);
    void                      test_binned_optimizer(void);
    void                      test_stacked_optimizer(void);
};

/***********************************************************************//**
 * @class TestGCTAPointing
 *
 * @brief Test suite for CTA pointing class testing
 *
 * This class defines a unit test suite for the CTA response class.
 ***************************************************************************/
class TestGCTAPointing : public GTestSuite {
public:

    // Constructors and destructors
    TestGCTAPointing(void) : GTestSuite() {}
    virtual ~TestGCTAPointing(void) {}

    // Methods
    virtual void              set(void);
    virtual TestGCTAPointing* clone(void) const;
    virtual std::string       classname(void) const { return "TestGCTAPointing"; }
    void                      test_load_table(void);
    void                      test_interpolate_altaz(void);
};

#endif /* TEST_CTA_HPP */
