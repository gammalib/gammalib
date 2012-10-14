/***************************************************************************
 *           GTestSuites.hpp  - Test Suites class for GammaLib             *
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
 * @file GTestSuites.hpp
 * @brief Test suites class definition
 * @author Jean-Baptiste Cayrou
 */

#ifndef GTESTSUITES_HPP
#define GTESTSUITES_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GBase.hpp"
#include "GLog.hpp"
#include "GXml.hpp"
#include "GTestSuite.hpp"


/***********************************************************************//**
 * @class GTestSuites
 *
 * @brief Interface for the test suites class
 *
 * This is a Test Suite container class.
 *
 * @todo Detailed description
 ***************************************************************************/
class GTestSuites : public GBase {
    
public:
    // Constructors and destructors
    GTestSuites(void);
    GTestSuites(const GTestSuites& suites);
    GTestSuites(const std::string& name);
    virtual ~GTestSuites(void);
    
    // Operators
    GTestSuites&      operator=(const GTestSuites& suites);
    GTestSuite&       operator[](const int& index);
    const GTestSuite& operator[](const int& index) const;
        
    // Methods
    void         clear(void);
    GTestSuites* clone(void) const;
    int          size(void) const;
    void         append(GTestSuite& suite);
    bool         run(void);
    void         save(const std::string& filename) const;
    std::string  name(void) const;
    void         name(const std::string& name);
    void         cout(bool cout);
    int          errors(void) const;
    int          failures(void) const;
    int          tests(void) const;
    time_t       timestamp(void) const;
    std::string  print(void) const;
    
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GTestSuites& suites);
    void free_members(void);
    void write(GXml& xml) const;
        
    // Protected members    
    std::string              m_name;       //!< Name of test suites
    std::vector<GTestSuite*> m_testsuites; //!< Vector of test suites
    time_t                   m_timestamp;  //!< Timestamp
    GLog                     m_log;        //!< Log
};

#endif /* GTESTSUITES_HPP */
