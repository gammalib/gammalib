/***************************************************************************
 *   GFitsTableCDoubleCol.i  - FITS table double precision complex column  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsTableCDoubleCol.i
 * @brief GFitsTableCDoubleCol class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsTableCDoubleCol.hpp"
%}


/***********************************************************************//**
 * @class GFitsTableCDoubleCol
 *
 * @brief SWIG interface for FITS table floating point column
 ***************************************************************************/
class GFitsTableCDoubleCol : public GFitsTableCol {
public:
    // Constructors and destructors
    GFitsTableCDoubleCol(void);
    GFitsTableCDoubleCol(const std::string& name, const int& length,
                         const int& size = 1);
    GFitsTableCDoubleCol(const GFitsTableCDoubleCol& column);
    virtual ~GFitsTableCDoubleCol(void);

    // Methods
    std::string     string(const int& row, const int& col = 0);
    double          real(const int& row, const int& col = 0);
    int             integer(const int& row, const int& col = 0);
    GFits::cdouble* data(void) { return m_data; }
    void            nulval(const GFits::cdouble* value);
    GFits::cdouble* nulval(void) { return m_nulval; }
};


/***********************************************************************//**
 * @brief GFitsTableCDoubleCol class extension
 ***************************************************************************/
%extend GFitsTableCDoubleCol {
    GFits::cdouble get(const int& row) {
        return (*self)(row);
    }
    GFits::cdouble get(const int& row, const int& col) {
        return (*self)(row, col);
    }
    void set(const int& row, const GFits::cdouble& value) {
        (*self)(row) = value;
    }
    void set(const int& row, const int& col, const GFits::cdouble& value) {
        (*self)(row, col) = value;
    }
    GFitsTableCDoubleCol copy() {
        return (*self);
    }
};
