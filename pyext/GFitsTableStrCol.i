/***************************************************************************
 *  GFitsTableStrCol.i  - FITS table string column class SWIG definition   *
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
 * @file GFitsTableStrCol.i
 * @brief GFitsTableStrCol class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsTableStrCol.hpp"
%}


/***********************************************************************//**
 * @class GFitsTableStrCol
 *
 * @brief SWIG interface for FITS table string column
 ***************************************************************************/
class GFitsTableStrCol : public GFitsTableCol {
public:
    // Constructors and destructors
    GFitsTableStrCol();
    GFitsTableStrCol(const std::string& name, const int& length,
                     const int& width, const int& size = 1);
    GFitsTableStrCol(const GFitsTableStrCol& column);
    ~GFitsTableStrCol();

    // Methods
    std::string  string(const int& row, const int& col = 0);
    double       real(const int& row, const int& col = 0);
    int          integer(const int& row, const int& col = 0);
    std::string* data(void);
    void         set_nullval(const std::string string);
};


/***********************************************************************//**
 * @brief GFitsTableStrCol class extension
 ***************************************************************************/
%extend GFitsTableStrCol {
    char *__str__() {
        static char str_buffer[1001];
        std::ostringstream buffer;
        buffer << *self;
	    std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 1001);
	    str_buffer[1000] = '\0';
	    return str_buffer;
    }
    std::string get(const int& row) {
        return (*self)(row);
    }
    std::string get(const int& row, const int& col) {
        return (*self)(row, col);
    }
    void set(const int& row, const std::string& value) {
        (*self)(row) = value;
    }
    void set(const int& row, const int& col, const std::string& value) {
        (*self)(row, col) = value;
    }
    GFitsTableStrCol copy() {
        return (*self);
    }
};
