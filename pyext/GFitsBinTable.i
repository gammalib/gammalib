/***************************************************************************
 *                GFitsBinTable.i - FITS binary table class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2012 by Juergen Knoedlseder                         *
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
 * @brief FITS binary table class definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsBinTable.hpp"
%}


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
    virtual void           clear(void);
    virtual GFitsBinTable* clone(void) const;
    HDUType                exttype(void) const { return HT_BIN_TABLE; }
};


/***********************************************************************//**
 * @brief GFitsBinTable class SWIG extension
 ***************************************************************************/
%extend GFitsBinTable {
    char *__str__() {
        return tochar(self->print());
    }
    GFitsBinTable copy() {
        return (*self);
    }
}
