/***************************************************************************
 *  GFitsTableShortCol.i  - FITS table short column class SWIG definition  *
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
 * @file GFitsTableShortCol.i
 * @brief GFitsTableShortCol class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsTableShortCol.hpp"
%}


/***********************************************************************//**
 * @class GFitsTableShortCol
 *
 * @brief SWIG interface for FITS table short integer column
 ***************************************************************************/
class GFitsTableShortCol : public GFitsTableCol {
public:
    // Constructors and destructors
    GFitsTableShortCol(void);
    GFitsTableShortCol(const std::string& name, const int& length,
                       const int& size = 1);
    GFitsTableShortCol(const GFitsTableShortCol& column);
    virtual ~GFitsTableShortCol(void);

    // Methods
    std::string string(const int& row, const int& col = 0);
    double      real(const int& row, const int& col = 0);
    int         integer(const int& row, const int& col = 0);
    short*      data(void) { return m_data; }
    void        nulval(const short* value);
    short*      nulval(void) { return m_nulval; }
};


/***********************************************************************//**
 * @brief GFitsTableShortCol class extension
 ***************************************************************************/
%extend GFitsTableShortCol {
    short __getitem__(int GFitsTableColInx[]) {
        if (GFitsTableColInx[0] == 1)
            return (*self)(GFitsTableColInx[1]);
        else if (GFitsTableColInx[0] == 2)
            return (*self)(GFitsTableColInx[1], GFitsTableColInx[2]);
        else {
            PyErr_SetString(PyExc_ValueError,"Too many arguments in tuple");
            return NULL;
        }
    }
    void __setitem__(int GFitsTableColInx[], short value) {
        if (GFitsTableColInx[0] == 1)
            (*self)(GFitsTableColInx[1]) = value;
        else if (GFitsTableColInx[0] == 2)
            (*self)(GFitsTableColInx[1], GFitsTableColInx[2]) = value;
        else {
            PyErr_SetString(PyExc_ValueError,"Too many arguments in tuple");
            return;
        }
    }
    GFitsTableShortCol copy() {
        return (*self);
    }
};
