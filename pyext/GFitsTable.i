/***************************************************************************
 *               GFitsTable.i - FITS abstract table base class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2013 by Juergen Knoedlseder                         *
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
 * @file GFitsTable.i
 * @brief FITS table abstract base class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsTable.hpp"
#include "GTools.hpp"
%}
/***********************************************************************//**
 * @brief Tuple to index conversion to provide column access.
 *
 * The following function provides conversion between a Python tuple and
 * an integer array. This allows table column access via tuples, such as in
 * a[(3,5,10)] = 10.0 or c = a[(2,9)]. Note that the typemap will be
 * globally defined after inclusing of this file, hence GFitsTable.i has
 * to be included before all table class swig files.
 ***************************************************************************/
%{
static int table_column_tuple(PyObject *input, int *ptr) {
    if (PySequence_Check(input)) {
        int size = PyObject_Length(input);
        if (size > 2) {
            PyErr_SetString(PyExc_ValueError,"Too many arguments in tuple");
            return 0;
        }
        ptr[0] = size;
        for (int i = 0; i < size; i++) {
            PyObject *o = PySequence_GetItem(input,i);
            if (!PyInt_Check(o)) {
                Py_XDECREF(o);
                PyErr_SetString(PyExc_ValueError,"Expecting a tuple of integers");
                return 0;
            }
            ptr[i+1] = (int)PyInt_AsLong(o);
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

// This is the typemap that makes use of the function defined above
%typemap(in) int GFitsTableColInx[ANY](int temp[3]) {
    if (!table_column_tuple($input,temp)) {
        return NULL;
    }
    $1 = &temp[0];
}

// This typecheck verifies that all arguments are integers. The typecheck
// is needed for using "int GFitsTableColInx" in overloaded methods.
%typemap(typecheck) int GFitsTableColInx[ANY] {
    $1 = 1;
    if (PySequence_Check($input)) {
        int size = PyObject_Length($input);
        for (int i = 0; i < size; i++) {
            PyObject *o = PySequence_GetItem($input,i);
            if (!PyInt_Check(o)) {
                $1 = 0;
                break;
            }
        }
    }
    else {
        if (!PyInt_Check($input)) {
            $1 = 0;
        }
    }
}


/***********************************************************************//**
 * @class GFitsTable
 *
 * @brief FITS table abstract base class Python interface definition
 ***************************************************************************/
class GFitsTable : public GFitsHDU {
public:
    // Constructors and destructors
    GFitsTable(void);
    explicit GFitsTable(const int& nrows);
    GFitsTable(const GFitsTable& table);
    virtual ~GFitsTable(void);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GFitsTable* clone(void) const = 0;
    virtual HDUType     exttype(void) const = 0;

    // Implemented Methods
    const int&     size(void) const;
    GFitsTableCol* set(const int& colnum, const GFitsTableCol& column);
    GFitsTableCol* set(const std::string& colname, const GFitsTableCol& column);
    GFitsTableCol* append(const GFitsTableCol& column);
    GFitsTableCol* insert(int colnum, const GFitsTableCol& column);
    GFitsTableCol* insert(const std::string& colname, const GFitsTableCol& column);
    void           append_rows(const int& nrows);
    void           insert_rows(const int& row, const int& nrows);
    void           remove_rows(const int& row, const int& nrows);
    const int&     nrows(void) const;
    const int&     ncols(void) const;
    bool           contains(const std::string& colname) const;
};


/***********************************************************************//**
 * @brief GFitsTable class extension
 ***************************************************************************/
%extend GFitsTable {
    GFitsTableCol* __getitem__(const int& colnum) {
        return (*self)[colnum];
    }
    GFitsTableCol* __getitem__(const std::string& colname) {
        return (*self)[colname];
    }
    void __setitem__(const int& colnum, const GFitsTableCol& col) {
        self->set(colnum, col);
        return;
    }
    void __setitem__(const std::string& colname, const GFitsTableCol& col) {
        self->set(colname, col);
        return;
    }
};
