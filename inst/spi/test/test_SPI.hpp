/***************************************************************************
 *                    test_SPI.hpp  -  Test SPI classes                    *
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
 * @file test_SPI.hpp
 * @brief Definition of SPI classes unit tests
 * @author Juergen Knoedlseder
 */

#ifndef TEST_SPI_HPP
#define TEST_SPI_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GSPILib.hpp"


/***********************************************************************//**
 * @class TestGSPIResponse
 *
 * @brief Test suite for SPI response class testing
 *
 * This class defines a unit test suite for the SPI response class.
 ***************************************************************************/
class TestGSPIResponse : public GTestSuite {
public:
    // Constructors and destructors
    TestGSPIResponse(void) : GTestSuite() {}
    virtual ~TestGSPIResponse(void) {}

    // Methods
    virtual void set(void);
    void         test_inst_dir(void);
    void         test_pointing(void);
    void         test_response(void);
};


/***********************************************************************//**
 * @class TestGSPIObservation
 *
 * @brief Test suite for SPI observation testing
 *
 * This class defines a unit test suite for testing of GSPIObservation.
 ***************************************************************************/
class TestGSPIObservation : public GTestSuite {
public:
    // Constructors and destructors
    TestGSPIObservation(void) : GTestSuite() {}
    virtual ~TestGSPIObservation(void) {}

    // Methods
    virtual void set(void);
    void         test_binned_obs(void);
    void         test_event_bin(void);
    void         test_event_cube(void);
};


/***********************************************************************//**
 * @class TestGSPIOptimize
 *
 * @brief Test suite for SPI optimizer testing
 *
 * This class defines a unit test suite for testing of SPI optimizers.
 ***************************************************************************/
class TestGSPIOptimize : public GTestSuite {
public:
    // Constructors and destructors
    TestGSPIOptimize(void) : GTestSuite() {}
    virtual ~TestGSPIOptimize(void) {}

    // Methods
    virtual void set(void);
    void         test_binned_optimizer(void);
};

#endif /* TEST_SPI_HPP */
