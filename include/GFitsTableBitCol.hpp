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
 * This class implements a FITS table Bit column.
 *
 * @todo Each Bit is actually stored in one Byte. To save memory a more
 * compact storage scheme could be implemented.
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
    char&             operator() (const int& row, const int& inx = 0);
    const char&       operator() (const int& row, const int& inx = 0) const;

    // Methods
    std::string string(const int& row, const int& col = 0);
    double      real(const int& row, const int& col = 0);
    int         integer(const int& row, const int& col = 0);
    char*       data(void) { return m_data; }
    void        nullval(const char* value);
    char*       nullval(void) { return m_nulval; }

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
    void              alloc_nulval(const char* value);
    void              init_data(void);
    void*             ptr_data(void) { return m_data; }
    void*             ptr_nulval(void) { return m_nulval; }

    // Private data area
    char* m_data;       //!< Data area
    char* m_nulval;     //!< NULL value
};

#endif /* GFITSTABLEBITCOL_HPP */
