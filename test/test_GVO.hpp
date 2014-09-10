/***************************************************************************
 *                   test_GVO.hpp - Test VO module                         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2014 by Juergen Knoedlseder                         *
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
 * @file test_GVO.hpp
 * @brief Definition of unit tests for VO module
 * @author Juergen Knoedlseder 
 */

#ifndef TEST_GVO_HPP
#define TEST_GVO_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"


/***********************************************************************//**
 * @class TestGVO
 *
 * @brief Test suite for VO module
 ***************************************************************************/
class TestGVO : public GTestSuite {

public:
    // Constructors and destructors
    TestGVO(void) : GTestSuite() {}
    virtual ~TestGVO(void) {}

    // Methods
    virtual void        set(void);
    virtual TestGVO*    clone(void) const;
    virtual std::string classname(void) const { return "TestGVO"; }
    void                test_GVOClient(void);

private:
    // Private members
};

#endif /* TEST_GVO_HPP */
