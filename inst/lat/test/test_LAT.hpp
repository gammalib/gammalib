/***************************************************************************
 *                    test_LAT.hpp  -  Test LAT classes                    *
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
 * @file test_LAT.hpp
 * @brief Definition of LAT classes unit tests
 * @author Juergen Knoedlseder
 */

#ifndef TEST_LAT_HPP
#define TEST_LAT_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"


/***********************************************************************//**
 * @class TestGLATResponse
 *
 * @brief Test suite for LAT response class testing
 *
 * This class defines a unit test suite for the LAT response class.
 ***************************************************************************/
class TestGLATResponse : public GTestSuite {
public:
    // Constructors and destructors
    TestGLATResponse(void) : GTestSuite() {}
    virtual ~TestGLATResponse(void) {}

    // Methods
    virtual void              set(void);
    virtual TestGLATResponse* clone(void) const;
    virtual std::string       classname(void) const { return "TestGLATResponse"; }
    void                      test_response_p6(void);
    void                      test_response_p7(void);
    void                      test_one_response(const std::string& irf);
};


/***********************************************************************//**
 * @class GLATLtCube
 *
 * @brief Test suite for LAT livetime cube
 ***************************************************************************/
class TestGLATLtCube : public GTestSuite {
public:
    // Constructors and destructors
    TestGLATLtCube(void) : GTestSuite() {}
    virtual ~TestGLATLtCube(void) {}

    // Methods
    virtual void            set(void);
    virtual TestGLATLtCube* clone(void) const;
    virtual std::string     classname(void) const { return "TestGLATLtCube"; }
    void                    test_ltcube_p6(void);
    void                    test_ltcube_p7(void);
    void                    test_one_ltcube(const std::string& datadir, const double& reference);
};


/***********************************************************************//**
 * @class TestGLATObservation
 *
 * @brief Test suite for LAT observation testing
 *
 * This class defines a unit test suite for testing of GLATObservation.
 ***************************************************************************/
class TestGLATObservation : public GTestSuite {
public:
    // Constructors and destructors
    TestGLATObservation(void) : GTestSuite() {}
    virtual ~TestGLATObservation(void) {}

    // Methods
    virtual void                 set(void);
    virtual TestGLATObservation* clone(void) const;
    virtual std::string          classname(void) const { return "TestGLATObservation"; }
    void                         test_unbinned_obs_p6(void);
    void                         test_unbinned_obs_p7(void);
    void                         test_binned_obs_p6(void);
    void                         test_binned_obs_p7(void);
    void                         test_one_unbinned_obs(const std::string& datadir);
    void                         test_one_binned_obs(const std::string& datadir, const std::string& irf);
};



/***********************************************************************//**
 * @class TestGLATOptimize
 *
 * @brief Test suite for LAT optimizer testing
 *
 * This class defines a unit test suite for testing of LAT optimizers.
 ***************************************************************************/
class TestGLATOptimize : public GTestSuite {
public:
    // Constructors and destructors
    TestGLATOptimize(void) : GTestSuite() {}
    virtual ~TestGLATOptimize(void) {}

    // Methods
    virtual void              set(void);
    virtual TestGLATOptimize* clone(void) const;
    virtual std::string       classname(void) const { return "TestGLATOptimize"; }
    void                      test_binned_optimizer_p6(void);
    void                      test_binned_optimizer_p7(void);
    void                      test_one_binned_optimizer(const std::string& datadir,
                                                        const std::string& irf,
                                                        const double*      fit_results);

};

#endif /* TEST_LAT_HPP */
