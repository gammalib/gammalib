/***************************************************************************
 *       GFitsBinTable.i  - FITS binary table class SWIG definition        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsBinTable.i
 * @brief Binary FITS table Python interface
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
    GFitsBinTable(void);
    GFitsBinTable(int nrows);
    GFitsBinTable(const GFitsBinTable& table);
    virtual ~GFitsBinTable(void);

    // Implemented pure virtual methods
    HDUType exttype(void) const;
};
