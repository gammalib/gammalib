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

#ifndef GFITSASCIITABLE_HPP
#define GFITSASCIITABLE_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsCfitsio.hpp"
#include "GFitsData.hpp"
#include "GFitsTableCol.hpp"

/* __ Namespaces _________________________________________________________ */


/***************************************************************************
 *                     GFitsAsciiTable class definition                    *
 ***************************************************************************/
class GFitsAsciiTable : public GFitsData {

    // I/O friends
    friend ostream& operator<< (ostream& os, const GFitsAsciiTable& table);

public:
    // Constructors and destructors
    GFitsAsciiTable();
    GFitsAsciiTable(const GFitsAsciiTable& table);
    virtual ~GFitsAsciiTable();

    // Operators
    GFitsAsciiTable& operator= (const GFitsAsciiTable& table);

    // Methods
    void             open(__fitsfile* fptr);
    void             close(void);
    GFitsAsciiTable* clone(void) const;
    GFitsTableCol*   column(const std::string colname);
    GFitsTableCol*   column(const int colnum);

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsAsciiTable& table);
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
GFitsAsciiTable* GFitsAsciiTable::clone(void) const 
{
    return new GFitsAsciiTable(*this);
}

#endif /* GFITSASCIITABLE_HPP */
