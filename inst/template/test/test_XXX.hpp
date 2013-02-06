/***************************************************************************
 *                    test_XXX.hpp  -  Test XXX classes                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
 * @file test_XXX.hpp
 * @brief Definition of XXX classes unit tests
 * @author Juergen Knoedlseder
 */

#ifndef TEST_XXX_HPP
#define TEST_XXX_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GXXXLib.hpp"


/***********************************************************************//**
 * @class TestGXXXResponse
 *
 * @brief Test suite for XXX response class testing
 *
 * This class defines a unit test suite for the XXX response class.
 ***************************************************************************/
class TestGXXXResponse : public GTestSuite {
public:
    // Constructors and destructors
    TestGXXXResponse(void) : GTestSuite() {}
    virtual ~TestGXXXResponse(void) {}

    // Methods
    virtual void set(void);
    void         test_inst_dir(void);
    void         test_pointing(void);
    void         test_response(void);
};


/***********************************************************************//**
 * @class TestGXXXObservation
 *
 * @brief Test suite for XXX observation testing
 *
 * This class defines a unit test suite for testing of GXXXObservation.
 ***************************************************************************/
class TestGXXXObservation : public GTestSuite {
public:
    // Constructors and destructors
    TestGXXXObservation(void) : GTestSuite() {}
    virtual ~TestGXXXObservation(void) {}

    // Methods
    virtual void set(void);
    void         test_binned_obs(void);
    void         test_event_bin(void);
    void         test_event_cube(void);
};


/***********************************************************************//**
 * @class TestGXXXOptimize
 *
 * @brief Test suite for XXX optimizer testing
 *
 * This class defines a unit test suite for testing of XXX optimizers.
 ***************************************************************************/
class TestGXXXOptimize : public GTestSuite {
public:
    // Constructors and destructors
    TestGXXXOptimize(void) : GTestSuite() {}
    virtual ~TestGXXXOptimize(void) {}

    // Methods
    virtual void set(void);
    void         test_binned_optimizer(void);
};

#endif /* TEST_XXX_HPP */
