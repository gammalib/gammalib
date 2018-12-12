/***************************************************************************
 *       GPythonOptimizerFunction.i - Python optimizer function class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2018 by Juergen Knoedlseder                              *
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
 * @file GPythonOptimizerFunction.i
 * @brief Python optimizer function class
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GOptimizerFunction.hpp"
#include "GTools.hpp"
#include "GVector.hpp"
#include "GMatrixSparse.hpp"
#include "GException.hpp"


/***********************************************************************//**
 * @class GPythonException
 *
 * @brief Python exception
 *
 * This class implements a Python exception. It allows to catch any error
 * message that occurs during function evaluation.
 ***************************************************************************/
class GPythonException : public GExceptionHandler {
public:
    class eval_error : public GExceptionHandler {
    public:
        eval_error(const std::string& origin) {
            
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
                char     *c_type   = NULL;
                char     *c_value  = NULL;
                if (PyUnicode_Check(py_type)) {
                    PyObject* temp_bytes = PyUnicode_AsEncodedString(py_type,
                                                                     "utf-8",
                                                                     "Error ~");
                    if (temp_bytes != NULL) {
                        c_type = PyString_AsString(temp_bytes);
                    }
                }
                else if (py_type != NULL) {
                    c_type = PyString_AsString(py_type);
                }
                if (c_type != NULL) {
                    m_message += std::string(c_type);
                    m_message += "\n";
                }
                if (PyUnicode_Check(py_value)) {
                    PyObject* temp_bytes = PyUnicode_AsEncodedString(py_value,
                                                                     "utf-8",
                                                                     "Error ~");
                    if (temp_bytes != NULL) {
                        c_value = PyString_AsString(temp_bytes);
                    }
                }
                else if (py_value != NULL) {
                    c_value = PyString_AsString(py_value);
                }
                if (c_value != NULL) {
                    m_message += std::string(c_value);
                }
                Py_DECREF(py_type);
                Py_DECREF(py_value);
            }
        }
    };
};


/***********************************************************************//**
 * @class GPythonOptimizerFunction
 *
 * @brief Python interface for optimizer function
 *
 * This class implements the call back for a Python optimizer function,
 * allowing to use the GammaLib optimizer classes from within Python.
 ***************************************************************************/
class GPythonOptimizerFunction : public GOptimizerFunction {
public:
    // Constructor
    GPythonOptimizerFunction(void) {
        m_py_object = NULL;
        m_value     = 0.0;
        m_gradient.clear();
        m_curvature.clear();
        m_pars.clear();
    }

    // Destructor
    virtual ~GPythonOptimizerFunction(void) {
        if (m_py_object != NULL) {
            Py_DECREF(m_py_object);
        }
    }

    // Evaluate Python function
    void eval(const GOptimizerPars& pars) {

        // Store parameters
        m_pars = pars;

        // Initialise value, gradient and curvature
        int npars = pars.size();
        if (npars < 1) {
            npars = 1;
        }

        // Initialise value, gradient vector and curvature matrix
        m_value     = 0.0;
        m_gradient  = GVector(npars);
        m_curvature = GMatrixSparse(npars,npars);

        // Execute evaluation
        PyObject* args = Py_BuildValue("()");
        PyObject* func = static_cast<PyObject*>(m_py_object);
        PyObject* res  = PyEval_CallObject(func, args);
        Py_DECREF(args);

        // Release result
        if (res != NULL) {
            Py_DECREF(res);
        }
        else {
            throw GPythonException::eval_error("GPythonOptimizerFunction::eval");
        }
    }

    // Return value
    double value(void) const {
        return m_value;
    }

    // Return gradient
    GVector* gradient(void) {
        return &m_gradient;
    }

    // Return curvature
    GMatrixSparse* curvature(void) {
        return &m_curvature;
    }

    // Return parameters
    const GOptimizerPars& _pars(void) const {
        return m_pars;
    }

    // Set Python callback function
    void _set_eval(PyObject* function) {
        m_py_object = function;
        Py_INCREF(function);
    }

    // Set value
    void _set_value(const double& value) {
        m_value = value;
    }

    // Callback information
    PyObject*      m_py_object;  //!< Python callback function
    double         m_value;      //!< Function value
    GVector        m_gradient;   //!< Function gradient
    GMatrixSparse  m_curvature;  //!< Function curvature
    GOptimizerPars m_pars;       //!< Current parameters
};
%}


/***********************************************************************//**
 * @class GPythonOptimizerFunction
 *
 * @brief Optimizer class Python interface definition
 ***************************************************************************/
class GPythonOptimizerFunction : public GOptimizerFunction {
public:        
    // Constructors and destructors
    GPythonOptimizerFunction(void);
    virtual ~GPythonOptimizerFunction(void);

    // Methods
    void           eval(const GOptimizerPars& pars);
    double         value(void) const;
    GVector*       gradient(void);
    GMatrixSparse* curvature(void);

    // Protected methods
    const GOptimizerPars& _pars(void) const;
    void                  _set_eval(PyObject* function);
    void                  _set_value(const double& value);
};
