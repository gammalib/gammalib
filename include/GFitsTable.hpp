/***************************************************************************
 *             GFitsTable.hpp - FITS table abstract base class             *
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
 * @file GFitsTable.hpp
 * @brief FITS table abstract base class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GFITSTABLE_HPP
#define GFITSTABLE_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsHDU.hpp"
#include "GFitsTableCol.hpp"


/***********************************************************************//**
 * @class GFitsTable
 *
 * @brief Abstract interface for FITS table
 *
 * This class defines the abstract interface for a FITS table. A FITS table
 * is a collection of columns with an identical number of rows. This class
 * provides high level access to table columns.
 *
 * @todo Implement remove_column method
 ***************************************************************************/
class GFitsTable : public GFitsHDU {

public:
    // Constructors and destructors
    GFitsTable(void);
    explicit GFitsTable(int nrows);
    GFitsTable(const GFitsTable& table);
    virtual ~GFitsTable(void);

    // Operators
    GFitsTable&          operator=(const GFitsTable& table);
    GFitsTableCol&       operator[](const int& colnum);
    const GFitsTableCol& operator[](const int& colnum) const;
    GFitsTableCol&       operator[](const std::string& colname);
    const GFitsTableCol& operator[](const std::string& colname) const;

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GFitsTable* clone(void) const = 0;
    virtual HDUType     exttype(void) const = 0;

    // Implemented Methods
    void        append_column(GFitsTableCol& column);
    void        insert_column(int colnum, GFitsTableCol& column);
    void        append_rows(const int& nrows);
    void        insert_rows(const int& rownum, const int& nrows);
    void        remove_rows(const int& rownum, const int& nrows);
    int         nrows(void) const;
    int         ncols(void) const;
    bool        hascolumn(const std::string& colname) const;
    std::string print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void  init_members(void);
    void  copy_members(const GFitsTable& table);
    void  free_members(void);
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
    GFitsTableCol*  alloc_column(int typecode) const;
    GFitsTableCol*  ptr_column(const std::string& colname) const;
};

#endif /* GFITSTABLE_HPP */
