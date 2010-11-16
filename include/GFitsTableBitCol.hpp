/***************************************************************************
 *            GFitsTableBitCol.hpp  - FITS table bit column class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsTableBitCol.hpp
 * @brief GFitsTableBitCol class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSTABLEBITCOL_HPP
#define GFITSTABLEBITCOL_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsTableCol.hpp"


/***********************************************************************//**
 * @class GFitsTableBitCol
 *
 * @brief Interface for FITS table Bit column
 *
 * This class implements a FITS table Bit column. Bits are stored internally
 * in an array of type char and is transferred to the file in junks of
 * 8 Bits (using the cfitsio type TBYTE).
 ***************************************************************************/
class GFitsTableBitCol : public GFitsTableCol {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GFitsTableBitCol& column);

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

    // Methods
    std::string    string(const int& row, const int& col = 0);
    double         real(const int& row, const int& col = 0);
    int            integer(const int& row, const int& col = 0);
    unsigned char* data(void) { return m_data; }
    void           nullval(const unsigned char* value);
    unsigned char* nullval(void) { return m_nulval; }

private:
    // Private methods
    void              init_members(void);
    void              copy_members(const GFitsTableBitCol& column);
    void              free_members(void);
    GFitsTableBitCol* clone(void) const;
    std::string       ascii_format(void) const;
    std::string       binary_format(void) const;
    void              alloc_data(void);
    void              release_data(void);
    void              alloc_nulval(const unsigned char* value);
    void              init_data(void);
    void*             ptr_data(void) { return m_data; }
    void*             ptr_nulval(void) { return m_nulval; }
    void              load_column(void);
    void              save_column(void);
    void              get_bit(const int& row, const int& inx);
    void              set_pending(void);

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
