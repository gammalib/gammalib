/***************************************************************************
 *     GFitsTable.i  - FITS table abstract base class SWIG interface       *
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
 * @file GFitsTable.i
 * @brief GFitsTable class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsTable.hpp"
%}


/***********************************************************************//**
 * @class GFitsTable
 *
 * @brief Abstract SWIG interface for the FITS table classes.
 ***************************************************************************/
class GFitsTable : public GFitsHDU {
public:
    // Constructors and destructors
    GFitsTable(void);
    GFitsTable(int nrows);
    GFitsTable(const GFitsTable& table);
    virtual ~GFitsTable(void);

    // Pure virtual methods
    virtual GFitsTable* clone(void) const = 0;

    // Implemented Methods
    virtual void           append_column(GFitsTableCol& column);
    virtual void           insert_column(int colnum, GFitsTableCol& column);
    virtual void           append_rows(const int& nrows);
    virtual void           insert_rows(const int& rownum, const int& nrows);
    virtual GFitsTableCol* column(const std::string& colname);
    virtual GFitsTableCol* column(const int& colnum);
    virtual int            nrows(void) const;
    virtual int            ncols(void) const;
};


/***********************************************************************//**
 * @brief GFitsTable class extension
 ***************************************************************************/
%extend GFitsTable {
    char *__str__() {
        static char str_buffer[100001];
        std::ostringstream buffer;
        buffer << *self;
	    std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 100001);
	    str_buffer[100000] = '\0';
	    return str_buffer;
    }
};
