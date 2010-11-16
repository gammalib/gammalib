/***************************************************************************
 *   GFitsTableLongCol.i  - FITS table long column class SWIG definiton    *
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
 * @file GFitsTableLongCol.i
 * @brief GFitsTableLongCol class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsTableLongCol.hpp"
%}


/***********************************************************************//**
 * @class GFitsTableLongCol
 *
 * @brief SWIG interface for FITS table long integer column
 ***************************************************************************/
class GFitsTableLongCol : public GFitsTableCol {
public:
    // Constructors and destructors
    GFitsTableLongCol(void);
    GFitsTableLongCol(const std::string& name, const int& length,
                      const int& size = 1);
    GFitsTableLongCol(const GFitsTableLongCol& column);
    virtual ~GFitsTableLongCol(void);

    // Methods
    std::string string(const int& row, const int& col = 0);
    double      real(const int& row, const int& col = 0);
    int         integer(const int& row, const int& col = 0);
    long*       data(void) { return m_data; }
    void        nulval(const long* value);
    long*       nulval(void) { return m_nulval; }
};


/***********************************************************************//**
 * @brief GFitsTableLongCol class extension
 ***************************************************************************/
%extend GFitsTableLongCol {
    char *__str__() {
        static char str_buffer[1001];
        std::ostringstream buffer;
        buffer << *self;
	    std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 1001);
	    str_buffer[1000] = '\0';
	    return str_buffer;
    }
    long get(const int& row) {
        return (*self)(row);
    }
    long get(const int& row, const int& col) {
        return (*self)(row, col);
    }
    void set(const int& row, const long& value) {
        (*self)(row) = value;
    }
    void set(const int& row, const int& col, const long& value) {
        (*self)(row, col) = value;
    }
    GFitsTableLongCol copy() {
        return (*self);
    }
};
