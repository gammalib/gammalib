/***************************************************************************
 *             GFitsTable.hpp  - FITS table abstract base class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2011 by Jurgen Knodlseder                           *
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
 * @file GFitsTable.hpp
 * @brief FITS table abstract base class interface definition
 * @author J. Knodlseder
 */

#ifndef GFITSTABLE_HPP
#define GFITSTABLE_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GLog.hpp"
#include "GFitsHDU.hpp"
#include "GFitsTableCol.hpp"


/* __ Structures _________________________________________________________ */


/***********************************************************************//**
 * @class GFitsTable
 *
 * @brief Abstract interface for FITS table
 *
 * This class defines the abstract interface for a FITS table. A FITS table
 * is a collection of columns with an identical number of rows.
 ***************************************************************************/
class GFitsTable : public GFitsHDU {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GFitsTable& table);
    friend GLog&         operator<< (GLog& log, const GFitsTable& table);

public:
    // Constructors and destructors
    GFitsTable(void);
    explicit GFitsTable(int nrows);
    GFitsTable(const GFitsTable& table);
    virtual ~GFitsTable(void);

    // Operators
    GFitsTable& operator= (const GFitsTable& table);

    // Pure virtual methods
    virtual GFitsTable* clone(void) const = 0;

    // Implemented Methods
    void           append_column(GFitsTableCol& column);
    void           insert_column(int colnum, GFitsTableCol& column);
    void           append_rows(const int& nrows);
    void           insert_rows(const int& rownum, const int& nrows);
    void           remove_rows(const int& rownum, const int& nrows);
    GFitsTableCol* column(const std::string& colname);
    GFitsTableCol* column(const int& colnum);
    int            nrows(void) const;
    int            ncols(void) const;
    std::string    print(void) const;

protected:
    // Protected methods
    void  data_open(void* vptr);
    void  data_save(void);
    void  data_close(void);
    void  data_connect(void* vptr);
    char* get_ttype(const int& colnum) const;
    char* get_tform(const int& colnum) const;
    char* get_tunit(const int& colnum) const;

    // Protected data area
    int             m_type;       //!< Table type (1=ASCII, 2=Binary)
    int             m_rows;       //!< Number of rows in table
    int             m_cols;       //!< Number of columns in table
    GFitsTableCol** m_columns;    //!< Array of table columns

private:
    // Private methods
    void            init_members(void);
    void            copy_members(const GFitsTable& table);
    void            free_members(void);
    GFitsTableCol*  alloc_column(int typecode) const;
};

#endif /* GFITSTABLE_HPP */
