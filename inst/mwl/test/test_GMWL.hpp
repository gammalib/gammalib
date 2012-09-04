/***************************************************************************
 *              test_MWL.cpp  -  Test multi-wavelength classes             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @file test_MWL.hpp
 * @brief Definition of multi-wavelength classes unit tests
 * @author Juergen Knoedlseder
 */

#ifndef TEST_GMWL_HPP
#define TEST_GMWL_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"


/***********************************************************************//**
 * @class TestGMWL
 *
 * @brief Test suite for multi-wavelength class testing
 *
 * This class defines a unit test suite for the multi-wavelength instrument
 * interface.
 ***************************************************************************/
class TestGMWL : public GTestSuite {
public:
    // Constructors and destructors
    TestGMWL(void) : GTestSuite() {}
    virtual ~TestGMWL(void) {}

    // Methods
    virtual void set(void);
    void         test_obs(void);
    void         test_optimizer(void);
};

#endif /* TEST_GMWL_HPP */
