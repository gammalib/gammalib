/***************************************************************************
 *                  test_COM.hpp  -  Test COMPTEL classes                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2017 by Juergen Knoedlseder                         *
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
 * @file test_COM.hpp
 * @brief Definition of COMTEL classes unit tests
 * @author Juergen Knoedlseder
 */

#ifndef TEST_COM_HPP
#define TEST_COM_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GCOMLib.hpp"
#include "GCOMSupport.hpp"


/***********************************************************************//**
 * @class TestGCOMSupport
 *
 * @brief Test suite for COMPTEL classes
 *
 * This class defines a unit test suite for the classes in the COMPTEL
 * instrument module.
 ***************************************************************************/
class TestGCOM : public GTestSuite {
public:
    // Constructors and destructors
    TestGCOM(void) : GTestSuite() {}
    virtual ~TestGCOM(void) {}

    // Methods
    virtual void        set(void);
    virtual TestGCOM*   clone(void) const;
    virtual std::string classname(void) const { return "TestGCOM"; }
    void                test_com_time(void);
    void                test_tim_class(void);
    void                test_oad_class(void);
    void                test_oads_class(void);
    void                test_inst_dir(void);
    void                test_response(void);
    void                test_binned_obs(void);
    void                test_event_bin(void);
    void                test_event_cube(void);
    void                test_binned_optimizer(void);
};

#endif /* TEST_COM_HPP */
