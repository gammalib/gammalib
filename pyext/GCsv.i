/***************************************************************************
 *               GCsv.i - Comma-separated values table class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2016 by Juergen Knoedlseder                         *
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
 * @file GCsv.hpp
 * @brief Comma-separated values table class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCsv.hpp"
#include "GTools.hpp"
%}
%include stl.i

/***********************************************************************//**
 * @brief Tuple to index conversion to provide element access.
 *
 * The following function provides conversion between a Python tuple and
 * an integer array. This allows table element access via tuples, such as in
 * a[3,5] = 10.0 or c = a[2,9]. Note that the typemap will be globally
 * defined after inclusing of this file.
 ***************************************************************************/
%{
static int csv_tuple(PyObject *input, int *ptr) {
    if (PySequence_Check(input)) {
        int size = PyObject_Length(input);
        if (size > 2) {
            PyErr_SetString(PyExc_ValueError,"Too many arguments in tuple");
            return 0;
        }
        for (int i = 0; i < size; i++) {
            PyObject *o = PySequence_GetItem(input,i);
            if (!PyInt_Check(o)) {
                Py_XDECREF(o);
                PyErr_SetString(PyExc_ValueError,"Expecting a tuple of integers");
                return 0;
            }
            ptr[i] = (int)PyInt_AsLong(o);
            Py_DECREF(o);
        }
        return 1;
    }
    else {
        ptr[0] = 1;
        if (!PyInt_Check(input)) {
            PyErr_SetString(PyExc_ValueError,"Expecting an integer");
            return 0;
        }
        ptr[1] = (int)PyInt_AsLong(input);
        return 1;       
    }
}
%}
%typemap(in) int GCsvInx[ANY](int temp[2]) {
   if (!csv_tuple($input,temp)) {
      return NULL;
   }
   $1 = &temp[0];
}


/***********************************************************************//**
 * @class GCsv
 *
 * @brief Comma-separated values table class definition
 ***************************************************************************/
class GCsv : public GBase {
public:
    // Constructors and destructors
    GCsv(void);
    GCsv(const int& nrows, const int& ncols);
    GCsv(const GFilename& filename, const std::string& sep = " ");
    GCsv(const GCsv& csv);
    virtual ~GCsv(void);
 
    // Methods
    void        clear(void);
    GCsv*       clone(void) const;
    std::string classname(void) const;
    int         size(void) const;
    const int&  ncols(void) const;
    const int&  nrows(void) const;
    const int&  precision(void) const;
    void        precision(const int& precision);
    std::string string(const int& row, const int& col) const;
    double      real(const int& row, const int& col) const;
    int         integer(const int& row, const int& col) const;
    void        string(const int& row, const int& col, const std::string& value);
    void        real(const int& row, const int& col, const double& value);
    void        integer(const int& row, const int& col, const int& value);
    void        load(const GFilename& filename, const std::string& sep = " ");
    void        save(const GFilename& filename, const std::string& sep = " ",
                     const bool& clobber = false) const;
};


/***********************************************************************//**
 * @brief GCsv class extension
 ***************************************************************************/
%extend GCsv {
    std::string __getitem__(int GCsvInx[]) {
        return (*self)(GCsvInx[0], GCsvInx[1]);
    }
    void __setitem__(int GCsvInx[], std::string value) {
        (*self)(GCsvInx[0], GCsvInx[1]) = value;
    }
    int __len__() {
        return (self->size());
    }
    GCsv copy() {
        return (*self);
    }
};
