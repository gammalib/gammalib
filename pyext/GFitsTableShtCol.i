/***************************************************************************
 *   GFitsTableShtCol.i  - FITS table short column class SWIG definition   *
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
 * @file GFitsTableShtCol.i
 * @brief GFitsTableShtCol class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsTableShtCol.hpp"
%}


/***********************************************************************//**
 * @class GFitsTableShtCol
 *
 * @brief SWIG interface for FITS table short integer column
 ***************************************************************************/
class GFitsTableShtCol : public GFitsTableCol {
public:
    // Constructors and destructors
    GFitsTableShtCol();
    GFitsTableShtCol(const std::string& name, const int& length,
                     const int& size = 1);
    GFitsTableShtCol(const GFitsTableShtCol& column);
    virtual ~GFitsTableShtCol();

    // Methods
    std::string string(const int& row, const int& col = 0);
    double      real(const int& row, const int& col = 0);
    int         integer(const int& row, const int& col = 0);
    short*      data(void);
    void        set_nullval(const short* value);
};


/***********************************************************************//**
 * @brief GFitsTableShtCol class extension
 ***************************************************************************/
%extend GFitsTableShtCol {
    char *__str__() {
        static char str_buffer[1001];
        std::ostringstream buffer;
        buffer << *self;
	    std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 1001);
	    str_buffer[1000] = '\0';
	    return str_buffer;
    }
    short get(const int& row) {
        return (*self)(row);
    }
    short get(const int& row, const int& col) {
        return (*self)(row, col);
    }
    void set(const int& row, const short& value) {
        (*self)(row) = value;
    }
    void set(const int& row, const int& col, const short& value) {
        (*self)(row, col) = value;
    }
    GFitsTableShtCol copy() {
        return (*self);
    }
};
