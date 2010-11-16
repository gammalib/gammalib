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
    char *__str__() {
        static char str_buffer[1001];
        std::ostringstream buffer;
        buffer << *self;
        std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 1001);
        str_buffer[1000] = '\0';
        return str_buffer;
    }
    double get(const int& row) {
        return (*self)(row);
    }
    double get(const int& row, const int& col) {
        return (*self)(row, col);
    }
    void set(const int& row, const double& value) {
        (*self)(row) = value;
    }
    void set(const int& row, const int& col, const double& value) {
        (*self)(row, col) = value;
    }
    GFitsTableDoubleCol copy() {
        return (*self);
    }
};
