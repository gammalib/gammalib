/***************************************************************************
 *                 test_SPI.hpp - Test INTEGRAL/SPI classes                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file test_SPI.cpp
 * @brief INTEGRAL/SPI test class definition
 * @author Juergen Knoedlseder
 */

#ifndef TEST_SPI_HPP
#define TEST_SPI_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GSPILib.hpp"


/***********************************************************************//**
 * @class TestGSPI
 *
 * @brief Test suite for testing of INTEGRAL/SPI classes
 *
 * This class defines a unit test suite for the INTEGRAL/SPI classes.
 ***************************************************************************/
class TestGSPI : public GTestSuite {
public:
    // Constructors and destructors
    TestGSPI(void) : GTestSuite() {}
    virtual ~TestGSPI(void) {}

    // Methods
    virtual void        set(void);
    virtual TestGSPI*   clone(void) const;
    virtual std::string classname(void) const { return "TestGSPI"; }
    void                test_obs(void);
    void                test_response(void);
};

#endif /* TEST_SPI_HPP */
