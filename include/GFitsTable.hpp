/***************************************************************************
 *             GFitsTable.hpp  - FITS table abstract base class            *
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
 * @file GFitsTable.hpp
 * @brief GFitsTable class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSTABLE_HPP
#define GFITSTABLE_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsCfitsio.hpp"
#include "GFitsData.hpp"
#include "GFitsTableCol.hpp"

/* __ Namespaces _________________________________________________________ */


/* __ Structures _________________________________________________________ */


/***********************************************************************//**
 * @class GFitsTable
 *
 * @brief Abstract interface for FITS table
 ***************************************************************************/
class GFitsTable : public GFitsData {

    // I/O friends
    friend ostream& operator<< (ostream& os, const GFitsTable& table);

public:
    // Constructors and destructors
    GFitsTable();
    GFitsTable(int nrows, int ncols);
    GFitsTable(const GFitsTable& table);
    virtual ~GFitsTable();

    // Operators
    GFitsTable& operator= (const GFitsTable& table);

    // Methods
    void           open(__fitsfile* fptr);
    void           save(void);
    void           close(void);
    GFitsTable*    clone(void) const = 0;
    GFitsTableCol* column(const std::string& colname);
    GFitsTableCol* column(const int& colnum);
    int            rows(void) const;
    int            cols(void) const;

protected:
    // Protected methods
    void  init_members(void);
    void  copy_members(const GFitsTable& table);
    void  free_members(void);
    void  connect(__fitsfile* fptr);
    char* get_ttype(const int& colnum) const;
    char* get_tform(const int& colnum) const;
    char* get_tunit(const int& colnum) const;

    // Protected data area
    int             m_type;       //!< Table type (1=ASCII, 2=Binary)
    int             m_rows;       //!< Number of rows in table
    int             m_cols;       //!< Number of columns in table
    GFitsTableCol** m_columns;    //!< Array of column pointers
};


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/

#endif /* GFITSTABLE_HPP */
