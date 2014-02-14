/***************************************************************************
 *                    test_CTA.hpp  -  Test CTA classes                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2014 by Juergen Knoedlseder                         *
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
    void                      test_response_aeff(void);
    void                      test_response_psf(void);
    void                      test_response_psf_king(void);
    void                      test_response_npsf(void);
    void                      test_response_edisp(void);
    void                      test_response_irf_diffuse(void);
    void                      test_response_npred_diffuse(void);
    void                      test_response(void);
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
    void                   test_modelbg_npred_xml(void);
    void                   test_modelbg_construct_fits(void);
    void                   test_model_inst_bgd(void);
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
    void                         test_unbinned_obs(void);
    void                         test_binned_obs(void);
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
    void                      test_unbinned_optimizer(void);
    void                      test_binned_optimizer(void);
};

#endif /* TEST_CTA_HPP */
