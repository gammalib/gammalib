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
/**
 * @file GFitsBinTable.hpp
 * @brief GFitsBinTable class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSBINTABLE_HPP
#define GFITSBINTABLE_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsCfitsio.hpp"
#include "GFitsData.hpp"
#include "GFitsTableCol.hpp"

/* __ Namespaces _________________________________________________________ */


/* __ Structures _________________________________________________________ */


/***********************************************************************//**
 * @class GFitsBinTable
 *
 * @brief Interface for FITS binary table
 ***************************************************************************/
class GFitsBinTable : public GFitsData {

    // I/O friends
    friend ostream& operator<< (ostream& os, const GFitsBinTable& table);

public:
    // Constructors and destructors
    GFitsBinTable();
    GFitsBinTable(int nrows, int ncols);
    GFitsBinTable(const GFitsBinTable& table);
    virtual ~GFitsBinTable();

    // Operators
    GFitsBinTable& operator= (const GFitsBinTable& table);

    // Methods
    void           open(__fitsfile* fptr);
    void           save(void);
    void           close(void);
    GFitsBinTable* clone(void) const;
    GFitsTableCol* column(const std::string& colname);
    GFitsTableCol* column(const int& colnum);

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsBinTable& table);
    void free_members(void);
    void connect(__fitsfile* fptr);

    // Private data area
    int             m_rows;       //!< Number of rows in table
    int             m_cols;       //!< Number of columns in table
    GFitsTableCol** m_columns;    //!< Array of column pointers
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
