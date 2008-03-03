/***************************************************************************
 *             GFitsAsciiTable.hpp  - FITS ASCII table class               *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsAsciiTable.hpp
 * @brief GFitsAsciiTable class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSASCIITABLE_HPP
#define GFITSASCIITABLE_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsCfitsio.hpp"
#include "GFitsData.hpp"
#include "GFitsTable.hpp"

/* __ Namespaces _________________________________________________________ */


/***********************************************************************//**
 * @class GFitsAsciiTable
 *
 * @brief Interface for FITS binary table
 ***************************************************************************/
class GFitsAsciiTable : public GFitsTable {

public:
    // Constructors and destructors
    GFitsAsciiTable();
    GFitsAsciiTable(int nrows, int ncols);
    GFitsAsciiTable(const GFitsAsciiTable& table);
    virtual ~GFitsAsciiTable();

    // Operators
    GFitsAsciiTable& operator= (const GFitsAsciiTable& table);

    // Methods
    GFitsAsciiTable* clone(void) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsAsciiTable& table);
    void free_members(void);

    // Private data area
};


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/

#endif /* GFITSASCIITABLE_HPP */
