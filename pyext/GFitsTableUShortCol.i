/***************************************************************************
 * GFitsTableUShortCol.i  - FITS table unsigned short column class SWIG def*
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
 * @file GFitsTableUShortCol.i
 * @brief GFitsTableUShortCol class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsTableUShortCol.hpp"
%}


/***********************************************************************//**
 * @class GFitsTableUShortCol
 *
 * @brief SWIG interface for FITS table short integer column
 ***************************************************************************/
class GFitsTableUShortCol : public GFitsTableCol {
public:
    // Constructors and destructors
    GFitsTableUShortCol(void);
    GFitsTableUShortCol(const std::string& name, const int& length,
                        const int& size = 1);
    GFitsTableUShortCol(const GFitsTableUShortCol& column);
    virtual ~GFitsTableUShortCol(void);

    // Methods
    std::string     string(const int& row, const int& col = 0);
    double          real(const int& row, const int& col = 0);
    int             integer(const int& row, const int& col = 0);
    unsigned short* data(void) { return m_data; }
    void            nulval(const unsigned short* value);
    unsigned short* nulval(void) { return m_nulval; }
};


/***********************************************************************//**
 * @brief GFitsTableUShortCol class extension
 ***************************************************************************/
%extend GFitsTableUShortCol {
    char *__str__() {
        static char str_buffer[1001];
        std::ostringstream buffer;
        buffer << *self;
	    std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 1001);
	    str_buffer[1000] = '\0';
	    return str_buffer;
    }
    unsigned short get(const int& row) {
        return (*self)(row);
    }
    unsigned short get(const int& row, const int& col) {
        return (*self)(row, col);
    }
    void set(const int& row, const unsigned short& value) {
        (*self)(row) = value;
    }
    void set(const int& row, const int& col, const unsigned short& value) {
        (*self)(row, col) = value;
    }
    GFitsTableUShortCol copy() {
        return (*self);
    }
};
