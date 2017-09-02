/***************************************************************************
 *                 test_XXX.hpp - Test [INSTRUMENT] classes                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) [YEAR] by [AUTHOR]                                       *
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
 * @file test_XXX.cpp
 * @brief [INSTRUMENT] test class definition
 * @author [AUTHOR]
 */

#ifndef TEST_XXX_HPP
#define TEST_XXX_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GXXXLib.hpp"


/***********************************************************************//**
 * @class TestGXXX
 *
 * @brief Test suite for testing of [INSTRUMENT] classes
 *
 * This class defines a unit test suite for the [INSTRUMENT] classes.
 ***************************************************************************/
class TestGXXX : public GTestSuite {
public:
    // Constructors and destructors
    TestGXXX(void) : GTestSuite() {}
    virtual ~TestGXXX(void) {}

    // Methods
    virtual void        set(void);
    virtual TestGXXX*   clone(void) const;
    virtual std::string classname(void) const { return "TestGXXX"; }
    void                test_obs(void);
    void                test_response(void);
};

#endif /* TEST_XXX_HPP */
