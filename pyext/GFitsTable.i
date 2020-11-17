/***************************************************************************
 *               GFitsTable.i - FITS abstract table base class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2020 by Juergen Knoedlseder                         *
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
 * @file GFitsTable.i
 * @brief FITS table abstract base class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsTable.hpp"
%}


/***********************************************************************//**
 * @class GFitsTable
 *
 * @brief FITS table abstract base class Python interface definition
 ***************************************************************************/
class GFitsTable : public GFitsHDU {
public:
    // Constructors and destructors
    GFitsTable(void);
    explicit GFitsTable(const int& nrows);
    GFitsTable(const GFitsTable& table);
    virtual ~GFitsTable(void);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GFitsTable* clone(void) const = 0;
    virtual std::string classname(void) const = 0;
    virtual HDUType     exttype(void) const = 0;

    // Implemented Methods
    GFitsTableCol* set(const int& colnum, const GFitsTableCol& column);
    GFitsTableCol* set(const std::string& colname, const GFitsTableCol& column);
    GFitsTableCol* append(const GFitsTableCol& column);
    GFitsTableCol* insert(int colnum, const GFitsTableCol& column);
    GFitsTableCol* insert(const std::string& colname, const GFitsTableCol& column);
    void           remove(const int& colnum);
    void           remove(const std::string& colname);
    void           append_rows(const int& nrows);
    void           insert_rows(const int& row, const int& nrows);
    void           remove_rows(const int& row, const int& nrows);
    const int&     nrows(void) const;
    const int&     ncols(void) const;
    bool           contains(const std::string& colname) const;
};


/***********************************************************************//**
 * @brief GFitsTable class extension
 ***************************************************************************/
%extend GFitsTable {
    GFitsTableCol* __getitem__(const int& colnum) {
        return (*self)[colnum];
    }
    GFitsTableCol* __getitem__(const std::string& colname) {
        return (*self)[colname];
    }
    void __setitem__(const int& colnum, const GFitsTableCol& col) {
        self->set(colnum, col);
        return;
    }
    void __setitem__(const std::string& colname, const GFitsTableCol& col) {
        self->set(colname, col);
        return;
    }
%pythoncode {
    def __getstate__(self):
        state = (gammalib.GFitsHDU.__getstate__(self), tuple([self[i] for i in range(self.ncols())]))
        return state
    def __setstate__(self, state):
        self.__init__()
        gammalib.GFitsHDU.__setstate__(self, state[0])
        for x in state[1]:
            self.append(x)
}
};
