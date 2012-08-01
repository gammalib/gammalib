/***************************************************************************
 *         GTestSuite.hpp  - Test Suite class for GammaLib                 *
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
 * @file GTestSuite.hpp
 * @brief Test Suite class definition
 * @author Jean-Baptiste Cayrou
 */

#ifndef GTESTSUITE_HPP
#define GTESTSUITE_HPP

/* __ Includes ___________________________________________________________ */

#include <iostream>
#include <vector>
#include <string>

#include "GLog.hpp"
#include "GException.hpp"
#include "GTestCase.hpp"

/* __ Forward declarations _______________________________________________ */
typedef void (*pfunction)(void);
class GTestSuite;

/***********************************************************************//**
 * @class GTestSuite
 *
 * @brief Test Suite class in order to perfom unit tests on GammaLib fixtures.
 *
 * TODO : explainations
 ***************************************************************************/
class GTestSuite{
    
    public:
        // Constructors and destructors
        GTestSuite(void);
        GTestSuite(const GTestSuite& testsuite);
        GTestSuite(const std::string& name);
        ~GTestSuite(void);

        // Operators
        GTestSuite& operator=(const GTestSuite& testsuite);
        GTestCase& operator[](const int& index);
        const GTestCase& operator[](const int& index) const;
        // Methods
        std::string name(void) const;
        void        name(const std::string& name);
        void        cout(bool cout);
        void        test_assert(bool result, const std::string& name);
        void        test_try(const std::string& name);
        void        test_try_success();
        void        test_try_failure(const std::string& message="",const std::string& type="");
        void        test_try_failure(const std::exception& e);
        void        test_function(pfunction function, const std::string& name);
        void        end_test();
        int         tests() const;
        int         errors() const;
        int         failures() const;
        int         success() const;
        std::time_t timestamp() const;
    
    // Protected methods
    protected:
        void        init_members(void);
        void        copy_members(const GTestSuite& testsuite);
        void        free_members(void);
    
    // Private members
    private:
        std::vector<GTestCase*> m_tests;    //!< Nomber of tests
        int                     m_failures; //!< Nomber of failures
        int                     m_errors;   //!< Nomber of errors
        std::string             m_name;     //!< Name of the test suite
        GLog                    m_log;      //!< Log
        std::time_t             m_timestamp;//!< Timestamp

};

#endif
