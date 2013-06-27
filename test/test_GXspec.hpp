/***************************************************************************
 *                 test_GXspec.hpp - Test Xspec module                     *
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
 * @file test_GXspec.hpp
 * @brief Definition of unit tests for Xspec module
 * @author Juergen Knoedlseder 
 */

#ifndef TEST_GXSPEC_HPP
#define TEST_GXSPEC_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"


/***********************************************************************//**
 * @class TestGXspec
 *
 * @brief Test suite for Xspec module
 ***************************************************************************/
class TestGXspec : public GTestSuite {

public:
    // Constructors and destructors
    TestGXspec(void) : GTestSuite() {}
    virtual ~TestGXspec(void) {}

    // Methods
    virtual void set(void);
    void         test_GPha(void);
    void         test_GArf(void);

private:
    // Private members
};

#endif /* TEST_GXSPEC_HPP */
