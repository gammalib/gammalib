/***************************************************************************
 *                    GTestCase.hpp - Test case class                      *
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
 * @file GTestCase.hpp
 * @brief Test case class interface definition
 * @author Jean-Baptiste Cayrou
 */
 
#ifndef GTESTCASE_H
#define GTESTCASE_H

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"


/***********************************************************************//**
 * @class GTestCase
 *
 * @brief Test class
 *
 * This class implements a single test case result. It contains information
 * about the test name, an optional message that is created in case of an
 * error or a failure, a type information (containing usually the name of
 * the class that is failing), and some information about the test duration.
 * Furthermore, the test case knows whether the test succeeded or failed,
 * and whether the test was an error or a failure test.
 ***************************************************************************/
class GTestCase : public GBase {
    
public:
    // Public enumerators
    enum ErrorKind {
        FAIL_TEST,
        ERROR_TEST
    };

    // Constructors and destructors
    GTestCase(void);
    GTestCase(const ErrorKind& kind, const std::string& name = "");
    GTestCase(const GTestCase& test);
    virtual ~GTestCase(void);

    // Operators
    GTestCase& operator=(const GTestCase& test);

    // Methods
    void               clear(void);
    GTestCase*         clone(void) const;
    std::string        classname(void) const;
    const std::string& name(void) const;
    void               name(const std::string& name);
    const std::string& message(void) const;
    void               message( const std::string& message);
    const std::string& type(void) const;
    void               type( const std::string& type);
    const ErrorKind&   kind(void) const;
    void               kind(const ErrorKind& kind);
    const bool&        has_passed(void) const;
    void               has_passed(const bool& has_passed);
    const double&      duration(void) const;
    void               duration(const double& duration);
    std::string        print(const GChatter& chatter = NORMAL) const;

protected:
     // Protected methods
    void init_members(void);
    void copy_members(const GTestCase& test);
    void free_members(void);
        
    // Private members
    std::string m_name;       //!< Test name
    std::string m_message;    //!< Test message
    std::string m_type;       //!< Test type
    bool        m_has_passed; //!< Boolean to check test success
    ErrorKind   m_kind;       //!< Kind of test case (failure or error test)
    double      m_duration;   //!< Test duration
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GTestCase").
 ***************************************************************************/
inline
std::string GTestCase::classname(void) const
{
    return ("GTestCase");
}


/***********************************************************************//**
 * @brief Return test case name
 *
 * @return Test case name.
 ***************************************************************************/
inline
const std::string& GTestCase::name(void) const
{
    return m_name;
}


/***********************************************************************//**
 * @brief Set test case name
 *
 * @param[in] name Test case name.
 ***************************************************************************/
inline
void GTestCase::name(const std::string& name)
{
    m_name = name;
    return;
}


/***********************************************************************//**
 * @brief Return test case message
 *
 * @return Test case message.
 ***************************************************************************/
inline
const std::string& GTestCase::message(void) const
{
    return m_message;
}


/***********************************************************************//**
 * @brief Set test case message
 *
 * @param[in] message Test case message.
 ***************************************************************************/
inline
void GTestCase::message(const std::string& message)
{
    m_message = message;
    return;
}


/***********************************************************************//**
 * @brief Return test case type
 *
 * @return Type of test case.
 ***************************************************************************/
inline
const std::string& GTestCase::type(void) const
{
    return m_type;
}


/***********************************************************************//**
 * @brief Set type of test case
 *
 * @param[in] type Type of test case.
 ***************************************************************************/
inline
void GTestCase::type(const std::string& type)
{
    m_type = type;
    return;
}


/***********************************************************************//**
 * @brief Return kind of test case
 *
 * @return Test case kind (FAIL_TEST or ERROR_TEST).
 *
 * Returns whether this test case is for failure testing (FAIL_TEST) or
 * error testing (ERROR_TEST).
 ***************************************************************************/
inline
const GTestCase::ErrorKind& GTestCase::kind(void) const
{
    return m_kind; 
}


/***********************************************************************//**
 * @brief Set kind of test case
 *
 * @param[in] kind Kind of test case (FAIL_TEST or ERROR_TEST)
 *
 * Specifies whether this test case is for failure testing (FAIL_TEST) or
 * error testing (ERROR_TEST).
 ***************************************************************************/
inline
void GTestCase::kind(const ErrorKind& kind)
{
    m_kind = kind;
    return;
}


/***********************************************************************//**
 * @brief Return whether the test passed
 *
 * @return True if test has passed, false otherwise.
 *
 * This method returns true if the test has passed, false otherwise.
 ***************************************************************************/
inline
const bool& GTestCase::has_passed(void) const
{
    return m_has_passed;
}


/***********************************************************************//**
 * @brief Set if test passed
 *
 * @param[in] has_passed Test passed (true or false)
 ***************************************************************************/
inline
void GTestCase::has_passed(const bool& has_passed)
{
    m_has_passed = has_passed;
    return;
}


/***********************************************************************//**
 * @brief Return test case duration
 *
 * @return Duration of test case.
 ***************************************************************************/
inline
const double& GTestCase::duration(void) const
{
    return m_duration;
}


/***********************************************************************//**
 * @brief Set test duration
 *
 * @param[in] duration Test duration.
 ***************************************************************************/
inline
void GTestCase::duration(const double& duration)
{
    m_duration = duration;
    return;
}

#endif /* GTESTCASE_H */
