/***************************************************************************
 *             test_GApplication.hpp - Test GApplication class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2013 by Jean-Baptiste Cayrou                        *
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
 * @file test_GApplication.hpp
 * @brief Definition of unit tests for GApplication classes
 * @author Jean-Baptiste Cayrou
 */

#ifndef TEST_GAPPLICATION_HPP
#define TEST_GAPPLICATION_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"


/***********************************************************************//**
 * @class TestGApplication
 *
 * @brief Test suite for GApplication class testing
 ***************************************************************************/
class TestGApplication : public GTestSuite {

public:
    // Constructors and destructors
    TestGApplication(void) : GTestSuite() {}
    virtual ~TestGApplication(void) {}

    // Methods
    virtual void              set(void);
    virtual TestGApplication* clone(void) const;
    virtual std::string       classname(void) const { return "TestGApplication"; }
    void                      test_constructor(void);
    void                      test_stream_logger(void);
    void                      test_C_logger(void);
    void                      test_GApplicationPar(void);
};

#endif /* TEST_GAPPLICATION_HPP */
