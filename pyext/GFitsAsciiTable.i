/***************************************************************************
 *      GFitsAsciiTable.hpp  - FITS ASCII table class SWIG definition      *
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
 * @file GFitsAsciiTable.i
 * @brief GFitsAsciiTable class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsAsciiTable.hpp"
%}
%feature("notabstract") GFitsAsciiTable;


/***********************************************************************//**
 * @class GFitsAsciiTable
 *
 * @brief SWIG interface for FITS binary table
 ***************************************************************************/
class GFitsAsciiTable : public GFitsTable {
public:
    // Constructors and destructors
    GFitsAsciiTable();
    GFitsAsciiTable(int nrows);
    GFitsAsciiTable(const GFitsAsciiTable& table);
    virtual ~GFitsAsciiTable();
};
