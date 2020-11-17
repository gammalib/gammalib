/***************************************************************************
 *      GFitsTableFloatCol.i - FITS table single precision column class    *
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
 * @file GFitsTableFloatCol.i
 * @brief FITS table single precision column class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsTableFloatCol.hpp"
%}


/***********************************************************************//**
 * @class GFitsTableFloatCol
 *
 * @brief FITS table single precision column class
 ***************************************************************************/
class GFitsTableFloatCol : public GFitsTableCol {
public:
    // Constructors and destructors
    GFitsTableFloatCol(void);
    GFitsTableFloatCol(const std::string& name, const int& nrows,
                       const int& size = 1);
    GFitsTableFloatCol(const GFitsTableFloatCol& column);
    virtual ~GFitsTableFloatCol(void);

    // Implement virtual methods
    virtual void                clear(void);
    virtual GFitsTableFloatCol* clone(void) const;
    virtual std::string         classname(void) const;
    virtual std::string         string(const int& row, const int& col = 0) const;
    virtual double              real(const int& row, const int& col = 0) const;
    virtual int                 integer(const int& row, const int& col = 0) const;
    virtual void                insert(const int& row, const int& nrows);
    virtual void                remove(const int& row, const int& nrows);
    virtual bool                is_loaded(void) const;
    
    // Other methods
    float* data(void);
    float* nulval(void);
    void   nulval(const float* value);
};


/***********************************************************************//**
 * @brief GFitsTableFloatCol class extension
 *
 * @todo We would like to return a reference in __getitem__ so that we can
 *       set the elements in an iterator. Declaring simply a reference
 *       shows an object pointer instead of the content. I still have to
 *       find out how to implement this (e.g. via typemaps).
 ***************************************************************************/
%extend GFitsTableFloatCol {
    float __getitem__(int GTuple1D2D[]) {
        if (GTuple1D2D[1] < 0 || GTuple1D2D[1] >= self->nrows()) {
            throw GException::out_of_range("__getitem__()", "Row index",
                                           GTuple1D2D[1], self->nrows());
        }
        if (GTuple1D2D[0] == 1) {
            return (*self)(GTuple1D2D[1]);
        }
        else {
            if (GTuple1D2D[2] >= 0 && GTuple1D2D[2] < self->elements(GTuple1D2D[1])) {
                return (*self)(GTuple1D2D[1], GTuple1D2D[2]);
            }
            else {
                throw GException::out_of_range("__getitem__()", "Column index",
                                               GTuple1D2D[2],
                                               self->elements(GTuple1D2D[1]));
            }
        }
    }
    void __setitem__(int GTuple1D2D[], float value) {
        if (GTuple1D2D[1] < 0 || GTuple1D2D[1] >= self->nrows()) {
            throw GException::out_of_range("__setitem__()", "Row index",
                                           GTuple1D2D[1], self->nrows());
        }
        if (GTuple1D2D[0] == 1) {
            (*self)(GTuple1D2D[1]) = value;
        }
        else {
            if (GTuple1D2D[2] >= 0 && GTuple1D2D[2] < self->elements(GTuple1D2D[1])) {
                (*self)(GTuple1D2D[1], GTuple1D2D[2]) = value;
            }
            else {
                throw GException::out_of_range("__setitem__()", "Column index",
                                               GTuple1D2D[2],
                                               self->elements(GTuple1D2D[1]));
            }
        }
    }
    GFitsTableFloatCol copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.name(), self.nrows(), self.number(), self.unit(),
                 [[self[row,inx] for inx in range(self.elements(row))] for row in range(self.nrows())])
        return state
    def __setstate__(self, state):
        self.__init__(state[0], state[1], state[2])
        self.unit(state[3])
        for row in range(len(state[4])):
            for inx in range(len(state[4][row])):
                self[row,inx] = state[4][row][inx]
}
};
