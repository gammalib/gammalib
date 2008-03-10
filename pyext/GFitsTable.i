/***************************************************************************
 *     GFitsTable.i  - FITS table abstract base class SWIG interface       *
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
class GFitsTable : public GFitsData {
public:
    // Constructors and destructors
    GFitsTable();
    GFitsTable(int nrows);
    GFitsTable(const GFitsTable& table);
    virtual ~GFitsTable();

    // Methods
    void           append_column(GFitsTableCol& column);
    void           insert_column(int colnum, GFitsTableCol& column);
    void           append_rows(const int& nrows);
    void           insert_rows(const int& rownum, const int& nrows);
    GFitsTableCol* column(const std::string& colname);
    GFitsTableCol* column(const int& colnum);
    int            nrows(void) const;
    int            ncols(void) const;

protected:
    GFitsTable* clone(void) const = 0;
};
