/***************************************************************************
 *            GTestSuite.i - Test suite class Python interface             *
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


/***********************************************************************//**
 * @class GTestSuite
 *
 * @brief Test suite Python interface defintion
 ***************************************************************************/
class GTestSuite {

public:
        
    // Constructors and destructors
    GTestSuite(void);
    GTestSuite(const GTestSuite& testsuite);
    GTestSuite(const std::string& name);
    virtual ~GTestSuite(void);

    // Methods
    void                      clear(void);
    int                       size(void) const;
    void                      append(pfunction function, const std::string& name);
    virtual void              set(void);
    bool                      run(void);
    std::string               name(void) const;
    void                      name(const std::string& name);
    void                      cout(bool cout);
    void                      test_assert(bool result,
                                          const std::string& name,
                                          const std::string& message = "");
    void                      test_value(double value,
                                         double value_expected,
                                         double eps = 0,
                                         const std::string& name = "",
                                         const std::string& message = "");
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
    void                      add_test(pfunction function, const std::string& name);
};


/***********************************************************************//**
 * @brief GTestSuite class extension
 ***************************************************************************/
%extend GTestSuite {
    char *__str__() {
        return tochar(self->print());
    }
    GTestCase& __getitem__(const int& index) {
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", index, self->size());
        }
    }
    void __setitem__(const int& index, const GTestCase& val) {
        if (index>=0 && index < self->size()) {
            (*self)[index] = val;
            return;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", index, self->size());
        }
    }
    void append(PyObject* function, const std::string& name) {
        // Push back the function pointer and the name
        self->m_functions.push_back(NULL);
        self->m_py_object.push_back(function);
        self->m_names.push_back(name);

        // Now register the function. We need this, otherwise we get
        // a Segmentation fault
        Py_INCREF(function);
    }
};
