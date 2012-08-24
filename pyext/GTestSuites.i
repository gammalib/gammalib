/***************************************************************************
 *   GTestSuites.i - Test suites class Python interface                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 Jean-Baptiste Cayrou                                *
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
 * @file GTestSuites.i
 * @brief Test suites class Python interface defintion
 * @author Jean-Baptiste Cayrou
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GTestSuites.hpp"
#include "GTools.hpp"
%}

// Include (int ARGC, char **ARGV) typemap to allow passing command line
// arguments to GApplication constructor
%include "argcargv.i"

/***********************************************************************//**
 * @class GTestSuites
 *
 * @brief Test suites Python interface defintion.
 ***************************************************************************/
class GTestSuites{
    
    public:
        // Constructors and destructors
        GTestSuites();
        GTestSuites(const GTestSuites& testsuites);
        GTestSuites(const std::string& name);
        ~GTestSuites();

        // Methods
        std::string name(void) const;
        void        name(const std::string& name);
        void        cout(bool cout);
        void        append(GTestSuite& testsuite);
        int         test_suites() const;
        int         errors(void) const;
        int         failures(void) const;
        int         tests(void) const;
        time_t      timestamp(void) const;
        bool        run(void);
        void        save(std::string filename);

};

/***********************************************************************//**
 * @brief GTestSuites class extension
 ***************************************************************************/
%extend GTestSuites {
    char *__str__() {
        return tochar(self->name());
    }

    GTestSuites& __getitem__(const int& index) {
        if (index >= 0 && index < self->test_suites()) {
            return (*self)[index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", index, self->test_suites());
        }
    }

    void __setitem__(const int& index, const GTestSuites& val) {
        if (index>=0 && index < self->test_suites()) {
            (*self)[index] = val;
            return;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", index, self->size());
        }
    }

    GTestSuites copy() {
        return (*self);
    }
};