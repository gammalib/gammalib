/***************************************************************************
 *        GFitsTableCol.hpp  - FITS table column abstract base class       *
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
 * @file GFitsTableCol.hpp
 * @brief FITS table column abstract base class definition
 * @author J. Knodlseder
 */

#ifndef GFITSTABLECOL_HPP
#define GFITSTABLECOL_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include <string>
#include "GLog.hpp"


/***********************************************************************//**
 * @class GFitsTableCol
 *
 * @brief Abstract interface for FITS table column
 *
 * This class implements a FITS table column. Vector columns are supported.
 ***************************************************************************/
class GFitsTableCol {

    // Friend classes
    friend class GFitsTable;

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GFitsTableCol& column);
    friend GLog&         operator<< (GLog& log, const GFitsTableCol& column);

public:
    // Constructors and destructors
    GFitsTableCol(void);
    explicit GFitsTableCol(const std::string& name, const int& length,
                           const int& number, const int& width);
    GFitsTableCol(const GFitsTableCol& column);
    virtual ~GFitsTableCol(void);

    // Operators
    GFitsTableCol& operator= (const GFitsTableCol& column);

    // Virtual Methods
    virtual std::string string(const int& row, const int& inx = 0) = 0;
    virtual double      real(const int& row, const int& inx = 0) = 0;
    virtual int         integer(const int& row, const int& inx = 0) = 0;
    virtual void        insert(const int& rownum, const int& nrows) = 0;
    virtual void        remove(const int& rownum, const int& nrows) = 0;

    // Base class Methods
    void        name(const std::string& name);
    void        unit(const std::string& unit);
    std::string name(void) const;
    std::string unit(void) const;
    int         colnum(void) const;
    int         type(void) const;
    int         repeat(void) const;
    int         width(void) const;
    int         number(void) const;
    int         length(void) const;
    int         anynul(void) const;
    std::string print(void) const;

protected:
    // Protected data area
    std::string m_name;      //!< Column name
    std::string m_unit;      //!< Column unit
    int         m_colnum;    //!< @brief Column number (starting from 1).
                             //!< This parameter is used to signal if a
                             //!< table column corresponds to a FITS file
                             //!< column. If it is set to 0 there is no
                             //!< correspondance.
    int         m_type;      //!< Column type
    int         m_repeat;    //!< Repeat value of column
    int         m_width;     //!< Width of single column element
    int         m_number;    //!< @brief Number of elements in column.
                             //!< m_number = m_repeat / m_width
    int         m_length;    //!< Length of column
    int         m_size;      //!< Size of allocated data area (0 if not loaded)
    int         m_anynul;    //!< Number of NULLs encountered
    void*       m_fitsfile;  //!< FITS file pointer associated with column

    // Protected virtual methods
    virtual GFitsTableCol* clone(void) const = 0;
    virtual std::string    ascii_format(void) const = 0;
    virtual std::string    binary_format(void) const = 0;
    virtual void           alloc_data(void) = 0;
    virtual void           init_data(void) = 0;
    virtual void*          ptr_data(void) = 0;
    virtual void*          ptr_nulval(void) = 0;

    // Protected methods
    virtual void           save(void);
    virtual void           fetch_data(void);
    virtual void           load_column(void);
    virtual void           save_column(void);
    virtual std::ostream&  dump_column(std::ostream& os) const;
    virtual GLog&          dump_column(GLog& log) const;
    virtual int            offset(const int& row, const int& inx) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsTableCol& column);
    void free_members(void);
    void connect(void* vptr);
};

#endif /* GFITSTABLECOL_HPP */
