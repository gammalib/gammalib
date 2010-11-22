/***************************************************************************
 * GFitsTableStringCol.i  - FITS table string column class SWIG definition *
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
 * @file GFitsTableStringCol.i
 * @brief GFitsTableStringCol class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsTableStringCol.hpp"
%}


/***********************************************************************//**
 * @class GFitsTableStringCol
 *
 * @brief SWIG interface for FITS table string column
 ***************************************************************************/
class GFitsTableStringCol : public GFitsTableCol {
public:
    // Constructors and destructors
    GFitsTableStringCol(void);
    GFitsTableStringCol(const std::string& name, const int& length,
                        const int& width, const int& size = 1);
    GFitsTableStringCol(const GFitsTableStringCol& column);
    virtual ~GFitsTableStringCol(void);

    // Methods
    std::string  string(const int& row, const int& col = 0);
    double       real(const int& row, const int& col = 0);
    int          integer(const int& row, const int& col = 0);
    std::string* data(void) { return m_data; }
    void         nulval(const std::string value);
    char*        nulval(void) { return m_nulval; }
};


/***********************************************************************//**
 * @brief GFitsTableStringCol class extension
 ***************************************************************************/
%extend GFitsTableStringCol {
    std::string __getitem__(int GFitsTableColInx[]) {
        if (GFitsTableColInx[0] == 1)
            return (*self)(GFitsTableColInx[1]);
        else
            return (*self)(GFitsTableColInx[1], GFitsTableColInx[2]);
    }
    void __setitem__(int GFitsTableColInx[], std::string value) {
        if (GFitsTableColInx[0] == 1)
            (*self)(GFitsTableColInx[1]) = value;
        else
            (*self)(GFitsTableColInx[1], GFitsTableColInx[2]) = value;
    }
    GFitsTableStringCol copy() {
        return (*self);
    }
};
