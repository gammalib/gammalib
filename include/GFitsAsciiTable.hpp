/***************************************************************************
 *             GFitsAsciiTable.hpp - FITS ASCII table class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2013 by Juergen Knoedlseder                         *
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
 * @file GFitsAsciiTable.hpp
 * @brief FITS ASCII table class definition
 * @author Juergen Knoedlseder
 */

#ifndef GFITSASCIITABLE_HPP
#define GFITSASCIITABLE_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsTable.hpp"


/***********************************************************************//**
 * @class GFitsAsciiTable
 *
 * @brief FITS ASCII table class
 *
 * The following ASCII table columns are supported: 
 * TSTRING (A),
 * TLONG (I),
 * TDOUBLE (F,D)
 * TFLOAT (E)
 ***************************************************************************/
class GFitsAsciiTable : public GFitsTable {

public:
    // Constructors and destructors
    GFitsAsciiTable(void);
    GFitsAsciiTable(const int& nrows);
    GFitsAsciiTable(const GFitsAsciiTable& table);
    virtual ~GFitsAsciiTable(void);

    // Operators
    GFitsAsciiTable& operator=(const GFitsAsciiTable& table);

    // Methods
    virtual void             clear(void);
    virtual GFitsAsciiTable* clone(void) const;
    HDUType                  exttype(void) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsAsciiTable& table);
    void free_members(void);
};


/***********************************************************************//**
 * @brief Return extension type
 *
 * @return Extension type (HT_ASCII_TABLE).
 ***************************************************************************/
inline
GFitsHDU::HDUType GFitsAsciiTable::exttype(void) const
{
    return (HT_ASCII_TABLE);
}

#endif /* GFITSASCIITABLE_HPP */
