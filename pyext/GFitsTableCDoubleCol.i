/***************************************************************************
 *   GFitsTableCDoubleCol.i  - FITS table double precision complex column  *
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
 * @file GFitsTableCDoubleCol.i
 * @brief FITS table double complex column class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsTableCDoubleCol.hpp"
%}


/***********************************************************************//**
 * @class GFitsTableCDoubleCol
 *
 * @brief FITS table double complex column Python interface definition
 ***************************************************************************/
class GFitsTableCDoubleCol : public GFitsTableCol {
public:
    // Constructors and destructors
    GFitsTableCDoubleCol(void);
    GFitsTableCDoubleCol(const std::string& name, const int& length,
                         const int& size = 1);
    GFitsTableCDoubleCol(const GFitsTableCDoubleCol& column);
    virtual ~GFitsTableCDoubleCol(void);

    // Implement virtual methods
    virtual void                  clear(void);
    virtual GFitsTableCDoubleCol* clone(void) const;
    virtual std::string           string(const int& row, const int& col = 0) const;
    virtual double                real(const int& row, const int& col = 0) const;
    virtual int                   integer(const int& row, const int& col = 0) const;
    virtual void                  insert(const int& rownum, const int& nrows);
    virtual void                  remove(const int& rownum, const int& nrows);
    
    // Other methods
    GFits::cdouble* data(void);
    void            nulval(const GFits::cdouble* value);
    GFits::cdouble* nulval(void);
};


/***********************************************************************//**
 * @brief GFitsTableCDoubleCol class extension
 *
 * @todo Implement __getitem__ correctly. This probably means that a real
 *       cdouble class is needed.
 *
 * @todo We would like to return a reference in __getitem__ so that we can
 *       set the elements in an iterator. Declaring simply a reference
 *       shows an object pointer instead of the content. I still have to
 *       find out how to implement this (e.g. via typemaps).
 ***************************************************************************/
%extend GFitsTableCDoubleCol {
/*
    GFits::cdouble __getitem__(int GFitsTableColInx[]) {
        if (GFitsTableColInx[0] == 1)
            return (*self)(GFitsTableColInx[1]);
        else
            return (*self)(GFitsTableColInx[1], GFitsTableColInx[2]);
    }
*/
    void __setitem__(int GFitsTableColInx[], GFits::cdouble value) {
        if (GFitsTableColInx[0] == 1)
            (*self)(GFitsTableColInx[1]) = value;
        else
            (*self)(GFitsTableColInx[1], GFitsTableColInx[2]) = value;
    }
    GFitsTableCDoubleCol copy() {
        return (*self);
    }
};
