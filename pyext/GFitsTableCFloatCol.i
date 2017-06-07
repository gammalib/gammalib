/***************************************************************************
 *   GFitsTableCFloatCol.i - FITS table single precision complex column    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2017 by Juergen Knoedlseder                         *
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
 * @file GFitsTableCFloatCol.i
 * @brief FITS table float complex column class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsTableCFloatCol.hpp"
%}


/***********************************************************************//**
 * @class GFitsTableCFloatCol
 *
 * @brief FITS table float complex column class
 ***************************************************************************/
class GFitsTableCFloatCol : public GFitsTableCol {
public:
    // Constructors and destructors
    GFitsTableCFloatCol(void);
    GFitsTableCFloatCol(const std::string& name, const int& nrows,
                        const int& size = 1);
    GFitsTableCFloatCol(const GFitsTableCFloatCol& column);
    virtual ~GFitsTableCFloatCol(void);

    // Implement virtual methods
    virtual void                 clear(void);
    virtual GFitsTableCFloatCol* clone(void) const;
    virtual std::string          classname(void) const;
    virtual std::string          string(const int& row, const int& col = 0) const;
    virtual double               real(const int& row, const int& col = 0) const;
    virtual int                  integer(const int& row, const int& col = 0) const;
    virtual void                 insert(const int& row, const int& nrows);
    virtual void                 remove(const int& row, const int& nrows);
    virtual bool                 is_loaded(void) const;
    
    // Other methods
    GFits::cfloat* data(void);
    GFits::cfloat* nulval(void);
    void           nulval(const GFits::cfloat* value);
};


/***********************************************************************//**
 * @brief GFitsTableCFloatCol class extension
 *
 * @todo Implement __getitem__ method. This probably means that a real
 *       cfloat class is needed.
 *
 * @todo We would like to return a reference in __getitem__ so that we can
 *       set the elements in an iterator. Declaring simply a reference
 *       shows an object pointer instead of the content. I still have to
 *       find out how to implement this (e.g. via typemaps).
 ***************************************************************************/
%extend GFitsTableCFloatCol {
/*
    GFits::cfloat __getitem__(int GFitsTableColInx[]) {
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
*/
    void __setitem__(int GFitsTableColInx[], GFits::cfloat value) {
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
    GFitsTableCFloatCol copy() {
        return (*self);
    }
};
