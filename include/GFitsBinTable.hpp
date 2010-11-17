/***************************************************************************
 *              GFitsBinTable.hpp  - FITS binary table class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2010 by Jurgen Knodlseder                           *
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
#include "GFitsTable.hpp"


/***********************************************************************//**
 * @class GFitsBinTable
 *
 * @brief Interface for FITS binary table
 ***************************************************************************/
class GFitsBinTable : public GFitsTable {

public:
    // Constructors and destructors
    GFitsBinTable(void);
    GFitsBinTable(int nrows);
    GFitsBinTable(const GFitsBinTable& table);
    virtual ~GFitsBinTable(void);

    // Operators
    GFitsBinTable& operator= (const GFitsBinTable& table);

    // Implemented pure virtual methods
    HDUType exttype(void) const { return HT_BIN_TABLE; }

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
