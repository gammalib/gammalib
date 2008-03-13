/***************************************************************************
 *  GFitsTableLogCol.i  - FITS table logical column class SWIG definition  *
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
 * @file GFitsTableLogCol.i
 * @brief GFitsTableLogCol class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsTableLogCol.hpp"
%}


/***********************************************************************//**
 * @class GFitsTableLogCol
 *
 * @brief SWIG interface for FITS table logical column
 ***************************************************************************/
class GFitsTableLogCol : public GFitsTableCol {
public:
    // Constructors and destructors
    GFitsTableLogCol();
    GFitsTableLogCol(const std::string& name, const int& length,
                     const int& size = 1);
    GFitsTableLogCol(const GFitsTableLogCol& column);
    virtual ~GFitsTableLogCol();

    // Methods
    std::string string(const int& row, const int& col = 0);
    double      real(const int& row, const int& col = 0);
    int         integer(const int& row, const int& col = 0);
    char*       data(void);
    void        set_nullval(const char* value);
};


/***********************************************************************//**
 * @brief GFitsTableLogCol class extension
 ***************************************************************************/
%extend GFitsTableLogCol {
    char *__str__() {
        static char str_buffer[1001];
        std::ostringstream buffer;
        buffer << *self;
	    std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 1001);
	    str_buffer[1000] = '\0';
	    return str_buffer;
    }
    char get(const int& row) {
        return (*self)(row);
    }
    char get(const int& row, const int& col) {
        return (*self)(row, col);
    }
    void set(const int& row, const char& value) {
        (*self)(row) = value;
    }
    void set(const int& row, const int& col, const char& value) {
        (*self)(row, col) = value;
    }
    GFitsTableLogCol copy() {
        return (*self);
    }
};
