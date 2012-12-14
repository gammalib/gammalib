/***************************************************************************
 *                test_GSupport.hpp - test support module                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Jean-Baptiste Cayrou                             *
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
 * @file test_GSupport.hpp
 * @brief Testing of support module
 * @author Jean-Baptiste Cayrou
 */

#ifndef TEST_GSUPPORT_HPP
#define TEST_GSUPPORT_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"



/***********************************************************************//**
 * @class TestGSupport
 *
 * @brief Test suite for support module
 ***************************************************************************/
class TestGSupport : public GTestSuite {

public:
    // Constructors and destructors
    TestGSupport(void) : GTestSuite(){ }
    virtual ~TestGSupport(void){ }

    // Methods
    virtual void set(void);
    void         test_expand_env(void);
    void         test_node_array(void);

private:
    // Private methods
    void test_node_array_interpolation(const int& num, const double* nodes);
};

#endif /* TEST_GSUPPORT_HPP */
