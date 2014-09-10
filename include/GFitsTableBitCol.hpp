/***************************************************************************
 *            GFitsTableBitCol.hpp - FITS table bit column class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2014 by Juergen Knoedlseder                         *
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
 * @author Juergen Knoedlseder
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
 * 8 Bits (using the CFITSIO type TBYTE).
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
    GFitsTableBitCol& operator=(const GFitsTableBitCol& column);
    bool&             operator()(const int& row, const int& inx = 0);
    const bool&       operator()(const int& row, const int& inx = 0) const;

    // Implemented virtual methods
    virtual void              clear(void);
    virtual GFitsTableBitCol* clone(void) const;
    virtual std::string       classname(void) const;
    virtual std::string       string(const int& row, const int& col = 0) const;
    virtual double            real(const int& row, const int& col = 0) const;
    virtual int               integer(const int& row, const int& col = 0) const;
    virtual void              insert(const int& row, const int& nrows);
    virtual void              remove(const int& row, const int& nrows);
    virtual bool              is_loaded(void) const;
    
    // Other methods
    unsigned char* data(void);
    unsigned char* nulval(void);
    void           nulval(const unsigned char* value);

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsTableBitCol& column);
    void free_members(void);
    void alloc_nulval(const unsigned char* value);
    void get_bit(const int& row, const int& inx);
    void set_pending(void);

    // Implemented virtual base class methods
    virtual void        alloc_data(void);
    virtual void        init_data(void);
    virtual void        fetch_data(void) const;
    virtual void        resize_data(const int& index, const int& number);
    virtual void        release_data(void);
    virtual void*       ptr_data(const int& index = 0);
    virtual void*       ptr_nulval(void);
    virtual std::string ascii_format(void) const;

    // Overloaded virtual methods
    virtual void load_column(void);
    virtual void save_column(void);

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


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GFitsTableBitCol").
 ***************************************************************************/
inline
std::string GFitsTableBitCol::classname(void) const
{
    return ("GFitsTableBitCol");
}


/***********************************************************************//**
 * @brief Checks if column has been loaded
 *
 * @return True if column has been loaded, false otherwise
 ***************************************************************************/
inline
bool GFitsTableBitCol::is_loaded(void) const
{
    return (m_data != NULL);
}


/***********************************************************************//**
 * @brief Returns pointer to column data
 *
 * @return Pointer to column data.
 ***************************************************************************/
inline
unsigned char* GFitsTableBitCol::data(void)
{
    return m_data;
}


/***********************************************************************//**
 * @brief Returns pointer to nul value
 *
 * @return Pointer to nul value.
 ***************************************************************************/
inline
unsigned char* GFitsTableBitCol::nulval(void)
{
    return m_nulval;
}


/***********************************************************************//**
 * @brief Returns void pointer to column data
 *
 * @return Void pointer to column data.
 ***************************************************************************/
inline
void* GFitsTableBitCol::ptr_data(const int& index)
{
    return (m_data+index);
}


/***********************************************************************//**
 * @brief Returns void pointer to nul value
 *
 * @return Void pointer to nul value.
 ***************************************************************************/
inline
void* GFitsTableBitCol::ptr_nulval(void)
{
    return m_nulval;
}

#endif /* GFITSTABLEBITCOL_HPP */
