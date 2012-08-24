/***************************************************************************
 *   GTestSuite.i - Test suite class Python interface                      *
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
 * @file GTestSuite.i
 * @brief Test suite class Python interface defintion
 * @author Jean-Baptiste Cayrou
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GTestSuite.hpp"
#include "GTools.hpp"
%}

// Include (int ARGC, char **ARGV) typemap to allow passing command line
// arguments to GApplication constructor
%include "argcargv.i"

/***********************************************************************//**
 * @class GTestSuite
 *
 * @brief Test suite Python interface defintion.
 ***************************************************************************/
class GTestSuite{

    public:
        // Constructors and destructors
                GTestSuite(void);
        GTestSuite(const GTestSuite& testsuite);
        GTestSuite(const std::string& name);
        virtual ~GTestSuite(void);

        // Methods
        virtual void                set(void) = 0;
        bool                        run(void);
        std::string                 name(void) const;
        void                        name(const std::string& name);
        void                        cout(bool cout);
        void                        test_assert(bool result,
                                                const std::string& name,
                                                const std::string& message="");
        void                        test_try(const std::string& name);
        void                        test_try_success();
        void                        test_try_failure(const std::string& message="",
                                                     const std::string& type="");
        void                        test_try_failure(const std::exception& e);
        GException::test_failure&   exception_failure(const std::string& message);
        GException::test_error&     exception_error(const std::string& message);
        void                        add_test(pfunction function, const std::string& name);
        int                         tests(void) const;
        int                         errors(void) const;
        int                         failures(void) const;
        int                         success(void) const;
        time_t                      timestamp(void) const;
        double                      duration(void) const;

};

/***********************************************************************//**
 * @brief GTestSuite class extension
 ***************************************************************************/
%extend GTestSuite {
    char *__str__() {
        return tochar(self->name());
    }

    GTestSuite& __getitem__(const int& index) {
        if (index >= 0 && index < self->tests()) {
            return (*self)[index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", index, self->tests());
        }
    }

    void __setitem__(const int& index, const GTestSuite& val) {
        if (index>=0 && index < self->test_suites()) {
            (*self)[index] = val;
            return;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", index, self->tests());
        }
    }

    GTestSuite copy() {
        return (*self);
    }

};