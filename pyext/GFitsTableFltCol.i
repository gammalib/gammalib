/***************************************************************************
 *   GFitsTableFltCol.i  - FITS table float column class SWIG definition   *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsTableFltCol.i
 * @brief GFitsTableFltCol class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsTableFltCol.hpp"
%}


/***********************************************************************//**
 * @class GFitsTableFltCol
 *
 * @brief SWIG interface for FITS table floating point column
 ***************************************************************************/
class GFitsTableFltCol : public GFitsTableCol {
public:
    // Constructors and destructors
    GFitsTableFltCol();
    GFitsTableFltCol(const std::string& name, const int& length,
                     const int& size = 1);
    GFitsTableFltCol(const GFitsTableFltCol& column);
    virtual ~GFitsTableFltCol();

    // Methods
    std::string string(const int& row, const int& col = 0);
    double      real(const int& row, const int& col = 0);
    int         integer(const int& row, const int& col = 0);
    float*      data(void);
    void        set_nullval(const float* value);
};


/***********************************************************************//**
 * @brief GFitsTableFltCol class extension
 ***************************************************************************/
%extend GFitsTableFltCol {
    char *__str__() {
        static char str_buffer[1001];
        std::ostringstream buffer;
        buffer << *self;
	    std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 1001);
	    str_buffer[1000] = '\0';
	    return str_buffer;
    }
    float get(const int& row) {
        return (*self)(row);
    }
    float get(const int& row, const int& col) {
        return (*self)(row, col);
    }
    void set(const int& row, const float& value) {
        (*self)(row) = value;
    }
    void set(const int& row, const int& col, const float& value) {
        (*self)(row, col) = value;
    }
    GFitsTableFltCol copy() {
        return (*self);
    }
};
