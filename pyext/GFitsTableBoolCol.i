/***************************************************************************
 * GFitsTableBoolCol.i  - FITS table boolean column class SWIG definition  *
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
 * @file GFitsTableBoolCol.i
 * @brief GFitsTableBoolCol class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsTableBoolCol.hpp"
%}


/***********************************************************************//**
 * @class GFitsTableBoolCol
 *
 * @brief SWIG interface for FITS table logical column
 ***************************************************************************/
class GFitsTableBoolCol : public GFitsTableCol {
public:
    // Constructors and destructors
    GFitsTableBoolCol(void);
    GFitsTableBoolCol(const std::string& name, const int& length,
                      const int& size = 1);
    GFitsTableBoolCol(const GFitsTableBoolCol& column);
    virtual ~GFitsTableBoolCol(void);

    // Methods
    std::string string(const int& row, const int& col = 0);
    double      real(const int& row, const int& col = 0);
    int         integer(const int& row, const int& col = 0);
    bool*       data(void) { return m_data; }
    void        nulval(const bool* value);
    bool*       nulval(void) { return m_nulval; }
};


/***********************************************************************//**
 * @brief GFitsTableBoolCol class extension
 ***************************************************************************/
%extend GFitsTableBoolCol {
    bool __getitem__(int GFitsTableColInx[]) {
        if (GFitsTableColInx[0] == 1)
            return (*self)(GFitsTableColInx[1]);
        else
            return (*self)(GFitsTableColInx[1], GFitsTableColInx[2]);
    }
    void __setitem__(int GFitsTableColInx[], bool value) {
        if (GFitsTableColInx[0] == 1)
            (*self)(GFitsTableColInx[1]) = value;
        else
            (*self)(GFitsTableColInx[1], GFitsTableColInx[2]) = value;
    }
    GFitsTableBoolCol copy() {
        return (*self);
    }
};
