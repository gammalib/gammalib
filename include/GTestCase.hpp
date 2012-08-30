/***************************************************************************
 *                           TestCase.hpp - Test Case Class                *
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
 * @file GTestCase.hpp
 * @brief Test case class definition
 * @author Jean-Baptiste Cayrou
 */
 
#ifndef GTESTCASE_H
#define GTESTCASE_H

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include <string>
#include <typeinfo>

#include "GTestSuite.hpp"

/* __ Forward declarations _______________________________________________ */
class GTestSuite;
typedef void (GTestSuite::*pfunction)(void);

/***********************************************************************//**
 * @class GTestCase
 *
 * @brief GTestCase class
 *
 * This class implements an test case object.
 ***************************************************************************/
class GTestCase{
    // Friend classes
    friend class GTestSuite;
    
    public:

        // public enumerators
        enum ErrorType{
            FAIL_TEST,
            ERROR_TEST
        };

        // Constructors and destructors
        GTestCase(void);
        GTestCase(const GTestCase& testcase);
        GTestCase(ErrorType type);
        GTestCase(const pfunction ptr_function, const std::string& name,GTestSuite* testsuite);
        GTestCase(ErrorType type, const std::string& name);
        ~GTestCase(void);
    
        // Operators
        GTestCase&  operator= (const GTestCase& testcase);
        
        // Methods
        std::string name(void) const;
        void        name(const std::string& name);
        std::string message(void) const;
        void        message( const std::string& message);
        std::string message_type(void) const;
        void        message_type( const std::string& message_type);
        pfunction   ptr_function(void) const;
        void        ptr_function(const pfunction);
        GTestSuite& testsuite(void) const;
        void        testsuite(GTestSuite* testsuite);
        ErrorType   type(void) const;
        void        type(ErrorType type);
        bool        is_passed() const;
        double      time() const;
        std::string print_result(void) const;
        void run(void);

     // Protected methods
    protected:
        void        init_members(void);
        void        copy_members(const GTestCase& testcase);
        void        free_members(void);
        
    // Private members
    private:
        std::string m_name; //!< name of the test
        std::string m_message; //!< message of the test
        std::string m_message_type; //!< type for the message
        pfunction   m_ptr_function; //!< function pointer
        GTestSuite* m_testsuite; //!< pointer to the test suite
        bool        m_passed; //!< boolean to check test success
        ErrorType   m_type; //!< type of the test case : failure or error test
        double      m_time; //!< duration of the test
};

#endif
