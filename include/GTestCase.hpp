/***************************************************************************
 *                    GTestCase.hpp  - Test case class                     *
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
#include "GBase.hpp"


/***********************************************************************//**
 * @class GTestCase
 *
 * @brief Test class
 *
 * This class implements a single test case result. It contains information
 * about the test name, an optional message that is created in case on an
 * error or a failure, a type information (containing usually the name of
 * the class that as failing), and some information about the test duration.
 * Furthermore, the test case knows whether the test succeeded or failed,
 * and whether the test was an error or a failure test.
 ***************************************************************************/
class GTestCase : public GBase {
    
public:
    // public enumerators
    enum ErrorKind {
        FAIL_TEST,
        ERROR_TEST
    };

    // Constructors and destructors
    GTestCase(void);
    GTestCase(const GTestCase& test);
    GTestCase(ErrorKind kind, const std::string& name = "");
    virtual ~GTestCase(void);

    // Operators
    GTestCase&  operator= (const GTestCase& test);

    // Methods
    void        clear(void);
    GTestCase*  clone(void) const;
    std::string name(void) const;
    void        name(const std::string& name);
    std::string message(void) const;
    void        message( const std::string& message);
    std::string type(void) const;
    void        type( const std::string& type);
    ErrorKind   kind(void) const;
    void        kind(ErrorKind kind);
    bool        passed(void) const;
    void        passed(const bool& passed);
    double      duration(void) const;
    void        duration(const double& duration);
    std::string print(void) const;

protected:
     // Protected methods
    void init_members(void);
    void copy_members(const GTestCase& test);
    void free_members(void);
        
    // Private members
    std::string m_name;        //!< Test name
    std::string m_message;     //!< Test message
    std::string m_type;        //!< Test type
    bool        m_passed;      //!< Boolean to check test success
    ErrorKind   m_kind;        //!< Kind of test case (failure or error test)
    double      m_duration;    //!< Test duration
};

#endif /* GTESTCASE_H */
