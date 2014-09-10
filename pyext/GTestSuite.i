/***************************************************************************
 *              GTestSuite.i - Abstract test suite base class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2014 Jean-Baptiste Cayrou                           *
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
 * @brief Abstract test suite base class definition
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
        m_py_names.clear();
        m_py_suites.clear();
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

    // Clone method (implementation makes this class non-abstract)
    GPythonTestSuite* clone(void) const {
        return new GPythonTestSuite(*this);
    }

    // Classname method (implementation makes this class non-abstract)
    std::string classname(void) const {
        return ("GPythonTestSuite");
    }

    // Generic test function for Python callback
    void test(void) {

        // Initialise failures, errors and test cases before testing
        int failures = m_py_suites[m_index]->failures();
        int errors   = m_py_suites[m_index]->errors();
        int cases    = m_py_suites[m_index]->size();

        // Execute unit test suite
        PyObject* args = Py_BuildValue("()");
        PyObject* func = static_cast<PyObject*>(m_py_objects[m_index]);
        PyObject* res  = PyEval_CallObject(func, args);
        Py_DECREF(args);

        // Handle any error that occured during the test
        if (res == NULL) {
            throw GPythonException::test_error("GPythonTestSuite::test");
        }
        else {
            Py_DECREF(res);
        }

        // Determine failures and errors after testing
        m_failures += (m_py_suites[m_index]->failures() - failures);
        m_errors   += (m_py_suites[m_index]->errors()   - errors);

        // Recover test cases
        for (int i = cases; i < m_py_suites[m_index]->size(); ++i) {
    
            // Get actual test case
            GTestCase&  test = (*m_py_suites[m_index])[i];

            // Kluge to set the test name correctly. For some reason, the
            // m_index points to the first test case when the Python
            // function is called, hence the name of the first test case
            // gets prepended. The name is thus replaced by the correct
            // name for all test cases following the first test case.
            if (m_index > 0) {
                std::string name = test.name();
                name.erase(0, m_names[0].length());
                name = m_names[m_index] + name;
                test.name(name);
            }

            // Clone the test case and push it on stack
            m_tests.push_back(new GTestCase(test));

        } // endfor: looped over all test cases

        // Return
        return;
    }

    // Append Python callback to function lists
    void append(PyObject* function, const std::string& name) {
        m_py_objects.push_back(function);
        m_py_names.push_back(name);
        m_py_suites.push_back(this);
        Py_INCREF(function);
        this->GTestSuite::append(static_cast<pfunction>(&GPythonTestSuite::test),
                                 name);
    }

    // Callback information
    std::vector<PyObject*>   m_py_objects; //!< Python callback function list
    std::vector<std::string> m_py_names;   //!< Python callback function names
    std::vector<GTestSuite*> m_py_suites;  //!< Python callback function pointers
};
%}


/***********************************************************************//**
 * @class GTestSuite
 *
 * @brief Abstract test suite Python interface definition
 ***************************************************************************/
class GTestSuite : public GBase {
public:        
    // Constructors and destructors
    GTestSuite(void);
    GTestSuite(const GTestSuite& testsuite);
    GTestSuite(const std::string& name);
    virtual ~GTestSuite(void);

    // Pure virtual methods
    virtual GTestSuite*       clone(void) const = 0;
    virtual void              set(void) = 0;

    // Methods
    void                      clear(void);
    int                       size(void) const;
    bool                      run(void);
    const std::string&        name(void) const;
    void                      name(const std::string& name);
    void                      cout(const bool& flag);
    void                      test_assert(const bool&        result,
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
    const int&                errors(void) const;
    const int&                failures(void) const;
    int                       success(void) const;
    const time_t&             timestamp(void) const;
    double                    duration(void) const;
};


/***********************************************************************//**
 * @brief GTestSuite class extension
 ***************************************************************************/
%extend GTestSuite {
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
 * @brief Test suite Python interface definition
 ***************************************************************************/
class GPythonTestSuite : public GTestSuite {
public:        
    // Constructors and destructors
    GPythonTestSuite(void);
    virtual ~GPythonTestSuite(void);

    // Methods
    void              set(void);
    GPythonTestSuite* clone(void) const;
    void              append(PyObject* function, const std::string& name);
};
