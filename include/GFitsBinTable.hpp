/***************************************************************************
 *              GFitsBinTable.hpp  - FITS binary table class               *
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

#ifndef GFITSBINTABLE_HPP
#define GFITSBINTABLE_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include "GFitsCfitsio.hpp"
#include "GFitsData.hpp"
#include "GFitsTableCol.hpp"

/* __ Namespaces _________________________________________________________ */


/* __ Structures _________________________________________________________ */


/***************************************************************************
 *                      GFitsBinTable class definition                     *
 ***************************************************************************/
class GFitsBinTable : public GFitsData {

public:
    // Constructors and destructors
    GFitsBinTable();
    GFitsBinTable(const GFitsBinTable& table);
    virtual ~GFitsBinTable();

    // Operators
    GFitsBinTable& operator= (const GFitsBinTable& table);

    // Methods
    void           open(__fitsfile*  fptr);
    void           close(void);
    GFitsBinTable* clone(void) const;
    GFitsTableCol* column(const std::string colname);
    GFitsTableCol* column(const int colnum);

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsBinTable& table);
    void free_members(void);

    // Private data area
    int             m_rows;
    int             m_cols;
    GFitsTableCol** m_columns;
};


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/
inline 
GFitsBinTable* GFitsBinTable::clone(void) const 
{
    return new GFitsBinTable(*this);
}

#endif /* GFITSBINTABLE_HPP */
