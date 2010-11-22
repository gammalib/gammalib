/***************************************************************************
 *  GFitsTableFloatCol.i  - FITS table float column class SWIG definition  *
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
 * @file GFitsTableFloatCol.i
 * @brief GFitsTableFloatCol class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsTableFloatCol.hpp"
%}


/***********************************************************************//**
 * @class GFitsTableFloatCol
 *
 * @brief SWIG interface for FITS table floating point column
 ***************************************************************************/
class GFitsTableFloatCol : public GFitsTableCol {
public:
    // Constructors and destructors
    GFitsTableFloatCol(void);
    GFitsTableFloatCol(const std::string& name, const int& length,
                       const int& size = 1);
    GFitsTableFloatCol(const GFitsTableFloatCol& column);
    virtual ~GFitsTableFloatCol(void);

    // Methods
    std::string string(const int& row, const int& col = 0);
    double      real(const int& row, const int& col = 0);
    int         integer(const int& row, const int& col = 0);
    float*      data(void) { return m_data; }
    void        nulval(const float* value);
    float*      nulval(void) { return m_nulval; }
};


/***********************************************************************//**
 * @brief GFitsTableFloatCol class extension
 ***************************************************************************/
%extend GFitsTableFloatCol {
    float __getitem__(int GFitsTableColInx[]) {
        if (GFitsTableColInx[0] == 1)
            return (*self)(GFitsTableColInx[1]);
        else
            return (*self)(GFitsTableColInx[1], GFitsTableColInx[2]);
    }
    void __setitem__(int GFitsTableColInx[], float value) {
        if (GFitsTableColInx[0] == 1)
            (*self)(GFitsTableColInx[1]) = value;
        else
            (*self)(GFitsTableColInx[1], GFitsTableColInx[2]) = value;
    }
    GFitsTableFloatCol copy() {
        return (*self);
    }
};
