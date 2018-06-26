/***************************************************************************
 *              GFitsAsciiTable.i - FITS ASCII table class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2018 by Juergen Knoedlseder                         *
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
 * @file GFitsAsciiTable.i
 * @brief FITS ASCII table class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsAsciiTable.hpp"
%}


/***********************************************************************//**
 * @class GFitsAsciiTable
 *
 * @brief SWIG interface for FITS ASCII table
 ***************************************************************************/
class GFitsAsciiTable : public GFitsTable {

public:
    // Constructors and destructors
    GFitsAsciiTable(void);
    explicit GFitsAsciiTable(const int& nrows);
    GFitsAsciiTable(const GFitsAsciiTable& table);
    virtual ~GFitsAsciiTable(void);

    // Methods
    virtual void             clear(void);
    virtual GFitsAsciiTable* clone(void) const;
    virtual std::string      classname(void) const;
    HDUType                  exttype(void) const;
};


/***********************************************************************//**
 * @brief GFitsAsciiTable class SWIG extension
 ***************************************************************************/
%extend GFitsAsciiTable {
    GFitsAsciiTable copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (gammalib.GFitsTable.__getstate__(self), self.nrows())
        return state
    def __setstate__(self, state):
        self.__init__(state[1])
        gammalib.GFitsTable.__setstate__(self, state[0])
}
}
