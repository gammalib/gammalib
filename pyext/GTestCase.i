/***************************************************************************
 *   GTestCase.i - Test case class Python interface                        *
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
 * @file GTestCase.i
 * @brief Test case class Python interface defintion
 * @author Jean-Baptiste Cayrou
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GTestCase.hpp"
#include "GTools.hpp"
%}

// Include (int ARGC, char **ARGV) typemap to allow passing command line
// arguments to GApplication constructor
%include "argcargv.i"

/***********************************************************************//**
 * @class GTestCase
 *
 * @brief Test case Python interface defintion.
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
};

/***********************************************************************//**
 * @brief GTestCase class extension
 ***************************************************************************/
%extend GTestCase {
    char *__str__() {
        return tochar(self->name());
    }

    GTestCase copy() {
        return (*self);
    }

};