/***************************************************************************
 *             GTestSuite.hpp - Abstract test suite base class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2014 by Jean-Baptiste Cayrou                        *
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
 * @brief Abstract test suite base class definition
 * @author Jean-Baptiste Cayrou
 */

#ifndef GTESTSUITE_HPP
#define GTESTSUITE_HPP

/* __ Forward declarations _______________________________________________ */
class GTestSuite;

/* __ Test suite pointer _________________________________________________ */
typedef void (GTestSuite::*pfunction)(void);

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GLog.hpp"
#include "GException.hpp"
#include "GTestCase.hpp"
#include "GTypemaps.hpp"


/***********************************************************************//**
 * @class GTestSuite
 *
 * @brief Abstract test suite class for unit testing on GammaLib fixtures
 *
 * @todo Detailed explanation.
 ***************************************************************************/
class GTestSuite : public GBase {
    
public:
    // Constructors and destructors
    GTestSuite(void);
    GTestSuite(const GTestSuite& testsuite);
    GTestSuite(const std::string& name);
    virtual ~GTestSuite(void);

    // Operators
    GTestSuite&      operator=(const GTestSuite& testsuite);
    GTestCase&       operator[](const int& index);
    const GTestCase& operator[](const int& index) const;

    // Pure virtual methods
    virtual GTestSuite*       clone(void) const = 0;
    virtual std::string       classname(void) const = 0;
    virtual void              set(void) = 0;

    // Other methods
    void                      clear(void);
    int                       size(void) const;
    void                      append(pfunction function, const std::string& name);
    bool                      run(void);
    const std::string&        name(void) const;
    void                      name(const std::string& name);
    void                      cout(const bool& flag);
    void                      test_assert(const bool&        result,
                                          const std::string& name,
                                          const std::string& message="");
    void                      test_value(const int&         value,
                                         const int&         expected,
                                         const std::string& name="",
                                         const std::string& message="");
    void                      test_value(const double&      value,
                                         const double&      expected,
                                         const double&      eps = 1.0e-10,
                                         const std::string& name="",
                                         const std::string& message="");
    void                      test_try(const std::string& name);
    void                      test_try_success(void);
    void                      test_try_failure(const std::string& message = "",
                                               const std::string& type = "");
    void                      test_try_failure(const std::exception& e);
    GException::test_failure& exception_failure(const std::string& message);
    GException::test_error&   exception_error(const std::string& message);
    const int&                errors(void) const;
    const int&                failures(void) const;
    int                       success(void) const;
    const time_t&             timestamp(void) const;
    double                    duration(void) const;
    std::string               print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const GTestSuite& testsuite);
    void        free_members(void);
    std::string format_name(const std::string& name);

    // Protected members
    std::string              m_name;       //!< Name of the test suite
    std::vector<std::string> m_names;      //!< Names of test functions
    std::vector<pfunction>   m_functions;  //!< Test functions of this suite
    std::vector<GTestCase*>  m_tests;      //!< List of test results
    std::vector<GTestCase*>  m_stack_try;  //!< Stack for nested try blocks
    int                      m_index;      //!< Index of actual test function
    int                      m_failures;   //!< Number of failures
    int                      m_errors;     //!< Number of errors
    GLog                     m_log;        //!< Log
    time_t                   m_timestamp;  //!< Timestamp
};


/***********************************************************************//**
 * @brief Return number of tests in test suite
 ***************************************************************************/
inline
int GTestSuite::size(void) const
{
    return m_tests.size();
}


/***********************************************************************//**
 * @brief Return test suite name
 ***************************************************************************/
inline
const std::string& GTestSuite::name(void) const
{
    return m_name;
}


/***********************************************************************//**
 * @brief Set Test Suite name
 *
 * @param[in] name Test suite name.
 ***************************************************************************/
inline
void GTestSuite::name(const std::string& name)
{
    m_name = name;
    return;
}


/***********************************************************************//**
 * @brief Enables/disables logging into standard output stream
 *
 * @param[in] flag Enable/disable logging (true/false).
 *
 * Enables or disables logging into the standard output stream.
 ***************************************************************************/
inline
void GTestSuite::cout(const bool& flag)
{
    m_log.cout(flag);
    return;
}


/***********************************************************************//**
 * @brief Return the number of errors
 *
 * @return Number of errors.
 ***************************************************************************/
inline
const int& GTestSuite::errors(void) const
{
    return m_errors; 
}


/***********************************************************************//**
 * @brief Return the number of failures
 *
 * @return Number of failures.
 ***************************************************************************/
inline
const int& GTestSuite::failures(void) const
{
    return m_failures; 
}


/***********************************************************************//**
 * @brief Return the timestamp
 *
 * @return Timestamp.
 *
 * The timestamp is set at the construction of the object.
 ***************************************************************************/
inline
const time_t& GTestSuite::timestamp(void) const
{
    return m_timestamp;
}

#endif /* GTESTSUITE_HPP */
