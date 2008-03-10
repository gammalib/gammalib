/***************************************************************************
 *    GFitsTableLngCol.i  - FITS table long column class SWIG definiton    *
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
 * @file GFitsTableLngCol.i
 * @brief GFitsTableLngCol class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsTableLngCol.hpp"
%}


/***********************************************************************//**
 * @class GFitsTableLngCol
 *
 * @brief SWIG interface for FITS table long integer column
 ***************************************************************************/
class GFitsTableLngCol : public GFitsTableCol {

    // I/O friends
    friend ostream& operator<< (ostream& os, const GFitsTableLngCol& column);

public:
    // Constructors and destructors
    GFitsTableLngCol();
    GFitsTableLngCol(const std::string& name, const int& length,
                     const int& size = 1);
    GFitsTableLngCol(const GFitsTableLngCol& column);
    virtual ~GFitsTableLngCol();

    // Methods
    std::string string(const int& row, const int& col = 0);
    double      real(const int& row, const int& col = 0);
    int         integer(const int& row, const int& col = 0);
    long*       data(void);
    void        set_nullval(const long* value);
};


/***********************************************************************//**
 * @brief GFitsTableLngCol class extension
 ***************************************************************************/
%extend GFitsTableLngCol {
    char *__str__() {
        static char str_buffer[100001];
        std::ostringstream buffer;
        buffer << *self;
	    std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 100001);
	    str_buffer[100000] = '\0';
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
    GFitsTableLngCol copy() {
        return (*self);
    }
};
