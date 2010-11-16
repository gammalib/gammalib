/***************************************************************************
 *  GFitsTableLongLongCol.i  - FITS table long long column class SWIG def. *
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
 * @file GFitsTableLongLongCol.i
 * @brief GFitsTableLongLongCol class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsTableLongLongCol.hpp"
%}


/***********************************************************************//**
 * @class GFitsTableLongLongCol
 *
 * @brief SWIG interface for FITS table long integer column
 ***************************************************************************/
class GFitsTableLongLongCol : public GFitsTableCol {
public:
    // Constructors and destructors
    GFitsTableLongLongCol(void);
    GFitsTableLongLongCol(const std::string& name, const int& length,
                          const int& size = 1);
    GFitsTableLongLongCol(const GFitsTableLongLongCol& column);
    virtual ~GFitsTableLongLongCol(void);

    // Methods
    std::string string(const int& row, const int& col = 0);
    double      real(const int& row, const int& col = 0);
    int         integer(const int& row, const int& col = 0);
    long long*  data(void) { return m_data; }
    void        nulval(const long long* value);
    long long*  nulval(void) { return m_nulval; }
};


/***********************************************************************//**
 * @brief GFitsTableLongLongCol class extension
 ***************************************************************************/
%extend GFitsTableLongLongCol {
    char *__str__() {
        static char str_buffer[1001];
        std::ostringstream buffer;
        buffer << *self;
	    std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 1001);
	    str_buffer[1000] = '\0';
	    return str_buffer;
    }
    long long get(const int& row) {
        return (*self)(row);
    }
    long long get(const int& row, const int& col) {
        return (*self)(row, col);
    }
    void set(const int& row, const long long& value) {
        (*self)(row) = value;
    }
    void set(const int& row, const int& col, const long long& value) {
        (*self)(row, col) = value;
    }
    GFitsTableLongLongCol copy() {
        return (*self);
    }
};
