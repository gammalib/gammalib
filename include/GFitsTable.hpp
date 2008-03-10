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
 *
 * This class implements a FITS table. A FITS table is a collection of
 * columns with an identical number of rows.
 ***************************************************************************/
class GFitsTable : public GFitsData {

    // I/O friends
    friend ostream& operator<< (ostream& os, const GFitsTable& table);

public:
    // Constructors and destructors
    GFitsTable();
    GFitsTable(int nrows);
    GFitsTable(const GFitsTable& table);
    virtual ~GFitsTable();

    // Operators
    GFitsTable& operator= (const GFitsTable& table);

    // Methods
    void           append_column(GFitsTableCol& column);
    void           insert_column(int colnum, GFitsTableCol& column);
    void           append_rows(const int& nrows);
    void           insert_rows(const int& rownum, const int& nrows);
    GFitsTableCol* column(const std::string& colname);
    GFitsTableCol* column(const int& colnum);
    int            nrows(void) const;
    int            ncols(void) const;

protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const GFitsTable& table);
    void        free_members(void);
    void        open(__fitsfile* fptr);
    void        save(void);
    void        close(void);
    void        connect(__fitsfile* fptr);
    GFitsTable* clone(void) const = 0;
    char*       get_ttype(const int& colnum) const;
    char*       get_tform(const int& colnum) const;
    char*       get_tunit(const int& colnum) const;

    // Protected data area
    int             m_type;       //!< Table type (1=ASCII, 2=Binary)
    int             m_rows;       //!< Number of rows in table
    int             m_cols;       //!< Number of columns in table
    GFitsTableCol** m_columns;    //!< Array of table columns
};


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/

#endif /* GFITSTABLE_HPP */
