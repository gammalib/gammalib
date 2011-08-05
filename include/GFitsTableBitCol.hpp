/***************************************************************************
 *            GFitsTableBitCol.hpp  - FITS table bit column class          *
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
 * @file GFitsTableBitCol.hpp
 * @brief FITS table bit column class interface definition
 * @author J. Knodlseder
 */

#ifndef GFITSTABLEBITCOL_HPP
#define GFITSTABLEBITCOL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GFitsTableCol.hpp"


/***********************************************************************//**
 * @class GFitsTableBitCol
 *
 * @brief FITS table Bit column
 *
 * This class implements a FITS table Bit column. Bits are stored internally
 * in an array of type char and is transferred to the file in junks of
 * 8 Bits (using the cfitsio type TBYTE).
 ***************************************************************************/
class GFitsTableBitCol : public GFitsTableCol {

public:
    // Constructors and destructors
    GFitsTableBitCol(void);
    GFitsTableBitCol(const std::string& name, const int& length,
                     const int& size = 1);
    GFitsTableBitCol(const GFitsTableBitCol& column);
    virtual ~GFitsTableBitCol(void);

    // Operators
    GFitsTableBitCol& operator= (const GFitsTableBitCol& column);
    bool&             operator() (const int& row, const int& inx = 0);
    const bool&       operator() (const int& row, const int& inx = 0) const;

    // Implement virtual methods
    virtual std::string string(const int& row, const int& col = 0) const;
    virtual double      real(const int& row, const int& col = 0) const;
    virtual int         integer(const int& row, const int& col = 0) const;
    virtual void        insert(const int& rownum, const int& nrows);
    virtual void        remove(const int& rownum, const int& nrows);
    
    // Other methods
    unsigned char* data(void) { return m_data; }
    void           nulval(const unsigned char* value);
    unsigned char* nulval(void) { return m_nulval; }

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsTableBitCol& column);
    void free_members(void);

    // Implemented pure virtual methods
    virtual GFitsTableBitCol* clone(void) const;
    virtual std::string       ascii_format(void) const;
    virtual std::string       binary_format(void) const;
    virtual void              alloc_data(void);
    virtual void              init_data(void);
    virtual void*             ptr_data(void) { return m_data; }
    virtual void*             ptr_nulval(void) { return m_nulval; }

    // Overloaded virtual methods
    virtual void load_column(void);
    virtual void save_column(void);

    // Other private methods
    void release_data(void);
    void alloc_nulval(const unsigned char* value);
    void get_bit(const int& row, const int& inx);
    void set_pending(void);

    // Private data area
    int            m_bits;           //!< Total number of Bits in column
    int            m_bytes_per_row;  //!< Number of Bytes per row
    int            m_bits_per_row;   //!< Number of Bits per row
    unsigned char* m_data;           //!< Data area
    unsigned char* m_nulval;         //!< NULL value

    // Bit access data area
    bool           m_bit_pending;    //!< Bit value has to be written back
    bool           m_bit_value;      //!< Actual bit to be accessed
    int            m_bit_byte;       //!< Row of actual bit to be accessed
    int            m_bit_mask;       //!< Index of actual bit to be accessed
};

#endif /* GFITSTABLEBITCOL_HPP */
