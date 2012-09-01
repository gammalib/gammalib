/***************************************************************************
 *           GTestSuite.hpp  - Test Suite class for GammaLib               *
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
 * @brief Test suite class definition
 * @author Jean-Baptiste Cayrou
 */

#ifndef GTESTSUITE_HPP
#define GTESTSUITE_HPP

/* __ Forward declarations _______________________________________________ */
class GTestSuite;

/* __ Test suite pointer _________________________________________________ */
typedef void (GTestSuite::*pfunction)(void);

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include <vector>
#include <string>
#include "GLog.hpp"
#include "GException.hpp"
#include "GTestCase.hpp"


/***********************************************************************//**
 * @class GTestSuite
 *
 * @brief Abstract test suite class for unit testing on GammaLib fixtures
 *
 * @todo Detailed explanation.
 ***************************************************************************/
class GTestSuite {
    
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
    virtual void set(void) = 0;

    // Other methods
    void                      clear(void);
    int                       size(void) const;
    void                      append(pfunction function, const std::string& name);
    virtual bool              run(void);
    std::string               name(void) const;
    void                      name(const std::string& name);
    void                      cout(bool cout);
    void                      test_assert(bool               result,
                                          const std::string& name,
                                          const std::string& message="");
    void                      test_value(const double&      value,
                                         const double&      expected,
                                         const double&      eps = 0.0,
                                         const std::string& name="",
                                         const std::string& message="");
    void                      test_try(const std::string& name);
    void                      test_try_success(void);
    void                      test_try_failure(const std::string& message = "",
                                               const std::string& type = "");
    void                      test_try_failure(const std::exception& e);
    GException::test_failure& exception_failure(const std::string& message);
    GException::test_error&   exception_error(const std::string& message);
    int                       errors(void) const;
    int                       failures(void) const;
    int                       success(void) const;
    time_t                    timestamp(void) const;
    double                    duration(void) const;

    // Old methods (will become obsolete)
    void add_test(pfunction function, const std::string& name);

protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const GTestSuite& testsuite);
    void        free_members(void);
    std::string format_name(const std::string& name);

    // Protected members
    std::string              m_name;       //!< Name of the test suite
    std::vector<pfunction>   m_functions;  //!< Test functions of this suite
    std::vector<std::string> m_names;      //!< Names of test functions
    std::vector<GTestCase*>  m_tests;      //!< List of test results
    std::vector<GTestCase*>  m_stack_try;  //!< Stack for nested try blocks
    int                      m_index;      //!< Index of actual test function
    int                      m_failures;   //!< Number of failures
    int                      m_errors;     //!< Number of errors
    GLog                     m_log;        //!< Log
    time_t                   m_timestamp;  //!< Timestamp

};

#endif /* GTESTSUITE_HPP */
