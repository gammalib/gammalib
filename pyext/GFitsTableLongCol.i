/***************************************************************************
 *       GFitsTableLongCol.i - FITS table long integer column class        *
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
 * @file GFitsTableLongCol.i
 * @brief FITS table long integer column class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsTableLongCol.hpp"
%}


/***********************************************************************//**
 * @class GFitsTableLongCol
 *
 * @brief FITS table long integer column class
 ***************************************************************************/
class GFitsTableLongCol : public GFitsTableCol {
public:
    // Constructors and destructors
    GFitsTableLongCol(void);
    GFitsTableLongCol(const std::string& name, const int& nrows,
                      const int& size = 1);
    GFitsTableLongCol(const GFitsTableLongCol& column);
    virtual ~GFitsTableLongCol(void);

    // Implement virtual methods
    virtual void               clear(void);
    virtual GFitsTableLongCol* clone(void) const;
    virtual std::string        classname(void) const;
    virtual std::string        string(const int& row, const int& col = 0) const;
    virtual double             real(const int& row, const int& col = 0) const;
    virtual int                integer(const int& row, const int& col = 0) const;
    virtual void               insert(const int& row, const int& nrows);
    virtual void               remove(const int& row, const int& nrows);
    virtual bool               is_loaded(void) const;
    
    // Other methods
    long* data(void);
    long* nulval(void);
    void  nulval(const long* value);
};


/***********************************************************************//**
 * @brief GFitsTableLongCol class extension
 *
 * @todo We would like to return a reference in __getitem__ so that we can
 *       set the elements in an iterator. Declaring simply a reference
 *       shows an object pointer instead of the content. I still have to
 *       find out how to implement this (e.g. via typemaps).
 ***************************************************************************/
%extend GFitsTableLongCol {
    long __getitem__(int GFitsTableColInx[]) {
        if (GFitsTableColInx[1] < 0 || GFitsTableColInx[1] >= self->nrows()) {
            throw GException::out_of_range("__getitem__()", "Row index",
                                           GFitsTableColInx[1], self->nrows());
        }
        if (GFitsTableColInx[0] == 1) {
            return (*self)(GFitsTableColInx[1]);
        }
        else {
            if (GFitsTableColInx[2] >= 0 &&
                GFitsTableColInx[2] < self->elements(GFitsTableColInx[1])) {
                return (*self)(GFitsTableColInx[1], GFitsTableColInx[2]);
            }
            else {
                throw GException::out_of_range("__getitem__()", "Column index",
                                               GFitsTableColInx[2],
                                               self->elements(GFitsTableColInx[1]));
            }
        }
    }
    void __setitem__(int GFitsTableColInx[], long value) {
        if (GFitsTableColInx[1] < 0 || GFitsTableColInx[1] >= self->nrows()) {
            throw GException::out_of_range("__setitem__()", "Row index",
                                           GFitsTableColInx[1], self->nrows());
        }
        if (GFitsTableColInx[0] == 1) {
            (*self)(GFitsTableColInx[1]) = value;
        }
        else {
            if (GFitsTableColInx[2] >= 0 && GFitsTableColInx[2] <
                self->elements(GFitsTableColInx[1])) {
                (*self)(GFitsTableColInx[1], GFitsTableColInx[2]) = value;
            }
            else {
                throw GException::out_of_range("__setitem__()", "Column index",
                                               GFitsTableColInx[2],
                                               self->elements(GFitsTableColInx[1]));
            }
        }
    }
    GFitsTableLongCol copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.name(), self.nrows(), self.number(),
                 [[self[row,inx] for inx in range(self.elements(row))] for row in range(self.nrows())])
        return state
    def __setstate__(self, state):
        self.__init__(state[0], state[1], state[2])
        for row in range(len(state[3])):
            for inx in range(len(state[3][row])):
                self[row,inx] = state[3][row][inx]
}
};
