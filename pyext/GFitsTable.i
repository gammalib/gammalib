/***************************************************************************
 *     GFitsTable.i  - FITS table abstract base class SWIG interface       *
 * ----------------------------------------------------------------------- *
 *  copyright : (C) 2008-2010 by Jurgen Knodlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsTable.i
 * @brief GFitsTable class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsTable.hpp"
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
%typemap(in) int GFitsTableColInx[ANY](int temp[3]) {
   if (!table_column_tuple($input,temp)) {
      return NULL;
   }
   $1 = &temp[0];
}


/***********************************************************************//**
 * @class GFitsTable
 *
 * @brief Abstract SWIG interface for the FITS table classes.
 ***************************************************************************/
class GFitsTable : public GFitsHDU {
public:
    // Constructors and destructors
    GFitsTable(void);
    explicit GFitsTable(int nrows);
    GFitsTable(const GFitsTable& table);
    virtual ~GFitsTable(void);

    // Pure virtual methods
    virtual GFitsTable* clone(void) const = 0;

    // Implemented Methods
    void           append_column(GFitsTableCol& column);
    void           insert_column(int colnum, GFitsTableCol& column);
    void           append_rows(const int& nrows);
    void           insert_rows(const int& rownum, const int& nrows);
    GFitsTableCol* column(const std::string& colname);
    GFitsTableCol* column(const int& colnum);
    int            nrows(void) const;
    int            ncols(void) const;
};


/***********************************************************************//**
 * @brief GFitsTable class extension
 ***************************************************************************/
%extend GFitsTable {
    char *__str__() {
        static std::string result = self->print();
        return ((char*)result.c_str());
    }
};
