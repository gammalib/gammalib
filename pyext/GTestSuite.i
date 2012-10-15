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
#include "GException.hpp"

/***********************************************************************//**
 * @class GPythonException
 *
 * @brief Exception for Python tests
 *
 * This class implements an exception for the Python tests. It allows to
 * catch any error message that occurs during Python unit testing.
 ***************************************************************************/
class GPythonException : public GExceptionHandler {
public:
    class test_error : public GExceptionHandler {
    public:
        test_error(const std::string& origin) {
            
            // Set origin and default message
            m_origin  = origin;
            m_message = "";

            // Set message
            if (PyErr_Occurred()) {

                // Fetch error type, value and traceback
                PyObject *type      = 0;
                PyObject *value     = 0;
                PyObject *traceback = 0;
                PyErr_Fetch(&type, &value, &traceback);

                // Clear error message
                PyErr_Clear();

                // Extract message strings
                PyObject *py_type  = PyObject_Str(type);
                PyObject *py_value = PyObject_Str(value);
                char     *c_type   = PyString_AsString(py_type);
                char     *c_value  = PyString_AsString(py_value);
                m_message += std::string(c_type);
                m_message += "\n";
                m_message += std::string(c_value);
                Py_DECREF(py_type);
                Py_DECREF(py_value);
            }
        }
    };
};


/***********************************************************************//**
 * @class GPythonTestSuite
 *
 * @brief Test suite for Python tests
 *
 * This is the test suite that should be used for Python tests. It derives
 * from the abstract GTestSuite base class and implements a generic test
 * method that is used as callback for Python testing. The class keeps an
 * array of PyObject pointers that is used to call the Python callbacks.
 * The generic test method profits from the fact that the function index is
 * a class member, hence the test method knows which Python callback function
 * it should execute.
 ***************************************************************************/
class GPythonTestSuite : public GTestSuite {
public:
    // Constructor
    GPythonTestSuite(void) {
        m_py_objects.clear();
    }

    // Destructor
    virtual ~GPythonTestSuite(void) {
        for (size_t i = 0; i < m_py_objects.size(); ++i) {
            Py_DECREF(m_py_objects[i]);
        }
    }

    // Dummy set method (implementation makes this class non-abstract)
    void set(void) {
    }

    // Generic test function for Python callback
    void test(void) {
        PyObject* args = Py_BuildValue("()");
        PyObject* func = static_cast<PyObject*>(m_py_objects[m_index]);
        PyObject* res  = PyEval_CallObject(func, args);
        Py_DECREF(args);
        if (res == NULL) {
            throw GPythonException::test_error("GPythonTestSuite::test");
        }
        else {
            Py_DECREF(res);
        }
        return;
    }

    // Append Python callback to function lists
    void append(PyObject* function, const std::string& name) {
        m_py_objects.push_back(function);
        Py_INCREF(function);
        add_test(static_cast<pfunction>(&GPythonTestSuite::test), name);
    }
    std::vector<PyObject*> m_py_objects; //!< Python callback function list
};
%}


/***********************************************************************//**
 * @class GTestSuite
 *
 * @brief Abstract test suite Python interface defintion
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
    virtual void              set(void) = 0;
    virtual bool              run(void);
    std::string               name(void) const;
    void                      name(const std::string& name);
    void                      cout(const bool& flag);
    void                      test_assert(bool result,
                                          const std::string& name,
                                          const std::string& message = "");
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
    int                       errors(void) const;
    int                       failures(void) const;
    int                       success(void) const;
    time_t                    timestamp(void) const;
    double                    duration(void) const;
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
};


/***********************************************************************//**
 * @class GTestSuite
 *
 * @brief Test suite Python interface defintion
 ***************************************************************************/
class GPythonTestSuite  : public GTestSuite {

public:
        
    // Constructors and destructors
    GPythonTestSuite(void);
    virtual ~GPythonTestSuite(void);

    // Methods
    void set(void);
    void append(PyObject* function, const std::string& name);
};
