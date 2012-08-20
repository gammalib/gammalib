/***************************************************************************
 *         GTestSuites.hpp  - Test Suites class for GammaLib               *
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
 * @brief Test Suites class definition
 * @author Jean-Baptiste Cayrou
 */

#ifndef GTESTSUITES_HPP
#define GTESTSUITES_HPP

#include <string>
#include <vector>

#include "GXml.hpp"
#include "GTools.hpp"
#include "GTestSuite.hpp"
#include "GTestCase.hpp"

/***********************************************************************//**
 * @class GTestSuites
 *
 * @brief Test Suites class 
 * This is a Test Suite container class.
 * TODO : Explainations
 ***************************************************************************/
class GTestSuites{
    
    public:
        // Constructors and destructors
        GTestSuites();
        GTestSuites(const GTestSuites& testsuites);
        GTestSuites(const std::string& name);
        ~GTestSuites();
    
        // Operators
        GTestSuites& operator=(const GTestSuites& testsuites);
        GTestSuite& operator[](const int& index);
        const GTestSuite& operator[](const int& index) const;
        
        // Methods
        std::string name() const;
        void        name(const std::string& name);
        void        cout(bool cout);
        void        append(GTestSuite& testsuite);
        int         test_suites() const;
        int         errors(void) const;
        int         failures(void) const;
        int         tests(void) const;
        std::time_t timestamp() const;
        bool        run();
        void        save(std::string filename);

    
    // Protected methods
    protected:
        void init_members();
        void copy_members(const GTestSuites& testsuites);
        void free_members();
        void write(GXml& xml);
        
    // Private members    
    private:
        std::string              m_name; //!< Test suites name
        std::vector<GTestSuite*> m_testsuites; //!< Test suite container
        std::time_t             m_timestamp;//!< Timestamp
        GLog                    m_log;      //!< Log
};

#endif
