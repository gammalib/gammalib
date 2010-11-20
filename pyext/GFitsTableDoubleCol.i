/***************************************************************************
 * GFitsTableDoubleCol.i  - FITS table double column class SWIG definiton  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsTableDoubleCol.i
 * @brief GFitsTableDoubleCol class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsTableDoubleCol.hpp"
%}


/***********************************************************************//**
 * @class GFitsTableDoubleCol
 *
 * @brief SWIG interface for FITS table double precision column
 ***************************************************************************/
class GFitsTableDoubleCol : public GFitsTableCol {
public:
    // Constructors and destructors
    GFitsTableDoubleCol(void);
    GFitsTableDoubleCol(const std::string& name, const int& length,
                        const int& size = 1);
    GFitsTableDoubleCol(const GFitsTableDoubleCol& column);
    virtual ~GFitsTableDoubleCol(void);

    // Methods
    std::string string(const int& row, const int& inx = 0);
    double      real(const int& row, const int& inx = 0);
    int         integer(const int& row, const int& inx = 0);
    double*     data(void) { return m_data; }
    void        nulval(const double* value);
    double*     nulval(void) { return m_nulval; }
};


/***********************************************************************//**
 * @brief GFitsTableDoubleCol class extension
 ***************************************************************************/
%extend GFitsTableDoubleCol {
   double __getitem__(int GFitsTableColInx[]) {
        if (GFitsTableColInx[0] == 1)
            return (*self)(GFitsTableColInx[1]);
        else if (GFitsTableColInx[0] == 2)
            return (*self)(GFitsTableColInx[1], GFitsTableColInx[2]);
        else {
            PyErr_SetString(PyExc_ValueError,"Too many arguments in tuple");
            return NULL;
        }
    }
    void __setitem__(int GFitsTableColInx[], double value) {
        if (GFitsTableColInx[0] == 1)
            (*self)(GFitsTableColInx[1]) = value;
        else if (GFitsTableColInx[0] == 2)
            (*self)(GFitsTableColInx[1], GFitsTableColInx[2]) = value;
        else {
            PyErr_SetString(PyExc_ValueError,"Too many arguments in tuple");
            return;
        }
    }
    GFitsTableDoubleCol copy() {
        return (*self);
    }
};
