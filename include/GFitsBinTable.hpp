/***************************************************************************
 *              GFitsBinTable.hpp  - FITS binary table class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsBinTable.hpp
 * @brief GFitsBinTable class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSBINTABLE_HPP
#define GFITSBINTABLE_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsData.hpp"
#include "GFitsTable.hpp"


/***********************************************************************//**
 * @class GFitsBinTable
 *
 * @brief Interface for FITS binary table
 ***************************************************************************/
class GFitsBinTable : public GFitsTable {

    // Friend classes
    friend class GFitsHDU;

public:
    // Constructors and destructors
    GFitsBinTable();
    GFitsBinTable(int nrows);
    GFitsBinTable(const GFitsBinTable& table);
    virtual ~GFitsBinTable();

    // Operators
    GFitsBinTable& operator= (const GFitsBinTable& table);

    // Methods

private:
    // Private methods
    void           init_members(void);
    void           copy_members(const GFitsBinTable& table);
    void           free_members(void);
    GFitsBinTable* clone(void) const;

    // Private data area
};

#endif /* GFITSBINTABLE_HPP */
