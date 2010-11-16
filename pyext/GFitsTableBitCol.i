/***************************************************************************
 *    GFitsTableBitCol.i  - FITS table bit column class SWIG definition    *
 * ----------------------------------------------------------------------- *
 *  copyright : (C) 2010 by Jurgen Knodlseder                              *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsTableBitCol.i
 * @brief GFitsTableBitCol class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsTableBitCol.hpp"
%}


/***********************************************************************//**
 * @class GFitsTableBitCol
 *
 * @brief SWIG interface for FITS table logical column
 ***************************************************************************/
class GFitsTableBitCol : public GFitsTableCol {
public:
    // Constructors and destructors
    GFitsTableBitCol(void);
    GFitsTableBitCol(const std::string& name, const int& length,
                     const int& size = 1);
    GFitsTableBitCol(const GFitsTableBitCol& column);
    virtual ~GFitsTableBitCol(void);

    // Methods
    std::string    string(const int& row, const int& col = 0);
    double         real(const int& row, const int& col = 0);
    int            integer(const int& row, const int& col = 0);
    unsigned char* data(void) { return m_data; }
    void           nullval(const unsigned char* value);
    unsigned char* nullval(void) { return m_nulval; }
};


/***********************************************************************//**
 * @brief GFitsTableBitCol class extension
 ***************************************************************************/
%extend GFitsTableBitCol {
    char *__str__() {
        static char str_buffer[1001];
        std::ostringstream buffer;
        buffer << *self;
	    std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 1001);
	    str_buffer[1000] = '\0';
	    return str_buffer;
    }
    bool get(const int& row) {
        return (*self)(row);
    }
    bool get(const int& row, const int& col) {
        return (*self)(row, col);
    }
    void set(const int& row, const bool& value) {
        (*self)(row) = value;
    }
    void set(const int& row, const int& col, const bool& value) {
        (*self)(row, col) = value;
    }
    GFitsTableBitCol copy() {
        return (*self);
    }
};
