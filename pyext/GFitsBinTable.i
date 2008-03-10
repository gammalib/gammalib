/***************************************************************************
 *       GFitsBinTable.i  - FITS binary table class SWIG definition        *
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
 * @file GFitsBinTable.i
 * @brief GFitsBinTable class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsBinTable.hpp"
%}
%feature("notabstract") GFitsBinTable;


/***********************************************************************//**
 * @class GFitsBinTable
 *
 * @brief SWIG interface for FITS binary table
 ***************************************************************************/
class GFitsBinTable : public GFitsTable {
public:
    // Constructors and destructors
    GFitsBinTable();
    GFitsBinTable(int nrows);
    GFitsBinTable(const GFitsBinTable& table);
    virtual ~GFitsBinTable();
};
