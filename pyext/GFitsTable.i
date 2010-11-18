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
    explicit GFitsTable(int nrows);
    GFitsTable(const GFitsTable& table);
    virtual ~GFitsTable(void);

    // Pure virtual methods
    virtual GFitsTable* clone(void) const = 0;

    // Implemented Methods
    void           append_column(GFitsTableCol& column);
    void           insert_column(int colnum, GFitsTableCol& column);
    void           append_rows(const int& nrows);
    void           insert_rows(const int& rownum, const int& nrows);
    GFitsTableCol* column(const std::string& colname);
    GFitsTableCol* column(const int& colnum);
    int            nrows(void) const;
    int            ncols(void) const;
};


/***********************************************************************//**
 * @brief GFitsTable class extension
 ***************************************************************************/
%extend GFitsTable {
    char *__str__() {
        std::string result = self->print();
        return ((char*)result.c_str());
    }
};
