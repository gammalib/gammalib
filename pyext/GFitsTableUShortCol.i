/***************************************************************************
 * GFitsTableUShortCol.i  - FITS table unsigned short column class SWIG def*
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
 * @file GFitsTableUShortCol.i
 * @brief GFitsTableUShortCol class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsTableUShortCol.hpp"
%}


/***********************************************************************//**
 * @class GFitsTableUShortCol
 *
 * @brief SWIG interface for FITS table short integer column
 ***************************************************************************/
class GFitsTableUShortCol : public GFitsTableCol {
public:
    // Constructors and destructors
    GFitsTableUShortCol(void);
    GFitsTableUShortCol(const std::string& name, const int& length,
                        const int& size = 1);
    GFitsTableUShortCol(const GFitsTableUShortCol& column);
    virtual ~GFitsTableUShortCol(void);

    // Methods
    std::string     string(const int& row, const int& col = 0);
    double          real(const int& row, const int& col = 0);
    int             integer(const int& row, const int& col = 0);
    unsigned short* data(void) { return m_data; }
    void            nulval(const unsigned short* value);
    unsigned short* nulval(void) { return m_nulval; }
};


/***********************************************************************//**
 * @brief GFitsTableUShortCol class extension
 ***************************************************************************/
%extend GFitsTableUShortCol {
    unsigned short __getitem__(int GFitsTableColInx[]) {
        if (GFitsTableColInx[0] == 1)
            return (*self)(GFitsTableColInx[1]);
        else if (GFitsTableColInx[0] == 2)
            return (*self)(GFitsTableColInx[1], GFitsTableColInx[2]);
        else {
            PyErr_SetString(PyExc_ValueError,"Too many arguments in tuple");
            return NULL;
        }
    }
    void __setitem__(int GFitsTableColInx[], unsigned short value) {
        if (GFitsTableColInx[0] == 1)
            (*self)(GFitsTableColInx[1]) = value;
        else if (GFitsTableColInx[0] == 2)
            (*self)(GFitsTableColInx[1], GFitsTableColInx[2]) = value;
        else {
            PyErr_SetString(PyExc_ValueError,"Too many arguments in tuple");
            return;
        }
    }
    GFitsTableUShortCol copy() {
        return (*self);
    }
};
