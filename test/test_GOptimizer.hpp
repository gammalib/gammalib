/***************************************************************************
 *               test_GOptimizer.hpp - test optimizer module               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2013 by Juergen Knoedlseder                         *
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
 * @file test_GOptimizer.hpp
 * @brief Definition of unit tests for optimizer module
 * @author Juergen Knoedlseder
 */

#ifndef TEST_GOPTIMIZER_HPP
#define TEST_GOPTIMIZER_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"


/***********************************************************************//**
 * @class TestGOptimizer
 *
 * @brief Test suite for optimizer module testing
 ***************************************************************************/
class TestGOptimizer : public GTestSuite {

public:
    // Constructors and destructors
    TestGOptimizer(void) : GTestSuite() {}
    virtual ~TestGOptimizer(void){}

    // Methods
    virtual void            set(void);
    virtual TestGOptimizer* clone(void) const;
    virtual std::string     classname(void) const { return "TestGOptimizer"; }
    void                    test_unbinned_optimizer(void);
    void                    test_binned_optimizer(void);
    void                    test_optimizer(const int& mode);
};

#endif /* TEST_GOPTIMIZER_HPP */
