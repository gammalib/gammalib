/***************************************************************************
 *               GFitsBinTable.hpp - FITS binary table class               *
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
 * @file GFitsBinTable.hpp
 * @brief FITS binary table class definition
 * @author Juergen Knoedlseder
 */

#ifndef GFITSBINTABLE_HPP
#define GFITSBINTABLE_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsTable.hpp"


/***********************************************************************//**
 * @class GFitsBinTable
 *
 * @brief FITS binary table class
 ***************************************************************************/
class GFitsBinTable : public GFitsTable {

public:
    // Constructors and destructors
    GFitsBinTable(void);
    GFitsBinTable(const int& nrows);
    GFitsBinTable(const GFitsBinTable& table);
    virtual ~GFitsBinTable(void);

    // Operators
    GFitsBinTable& operator=(const GFitsBinTable& table);

    // Implemented pure virtual methods
    virtual void           clear(void);
    virtual GFitsBinTable* clone(void) const;
    HDUType                exttype(void) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsBinTable& table);
    void free_members(void);
};


/***********************************************************************//**
 * @brief Return extension type
 *
 * @return Extension type (HT_BIN_TABLE).
 ***************************************************************************/
inline
GFitsHDU::HDUType GFitsBinTable::exttype(void) const
{
    return (HT_BIN_TABLE);
}

#endif /* GFITSBINTABLE_HPP */
