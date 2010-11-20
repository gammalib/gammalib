/***************************************************************************
 *         GFitsTableByteCol.hpp  - FITS table Byte column class           *
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
 * @file GFitsTableByteCol.hpp
 * @brief GFitsTableByteCol class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSTABLEBYTECOL_HPP
#define GFITSTABLEBYTECOL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GFitsTableCol.hpp"


/***********************************************************************//**
 * @class GFitsTableByteCol
 *
 * @brief Interface for FITS table Byte column
 *
 * This class implements a FITS table Byte column.
 ***************************************************************************/
class GFitsTableByteCol : public GFitsTableCol {

public:
    // Constructors and destructors
    GFitsTableByteCol(void);
    GFitsTableByteCol(const std::string& name, const int& length, const int& size = 1);
    GFitsTableByteCol(const GFitsTableByteCol& column);
    virtual ~GFitsTableByteCol(void);

    // Operators
    GFitsTableByteCol&   operator= (const GFitsTableByteCol& column);
    unsigned char&       operator() (const int& row, const int& inx = 0);
    const unsigned char& operator() (const int& row, const int& inx = 0) const;

    // Methods
    std::string    string(const int& row, const int& col = 0);
    double         real(const int& row, const int& col = 0);
    int            integer(const int& row, const int& col = 0);
    unsigned char* data(void) { return m_data; }
    void           nulval(const unsigned char* value);
    unsigned char* nulval(void) { return m_nulval; }

private:
    // Private methods
    void               init_members(void);
    void               copy_members(const GFitsTableByteCol& column);
    void               free_members(void);
    GFitsTableByteCol* clone(void) const;
    std::string        ascii_format(void) const;
    std::string        binary_format(void) const;
    void               alloc_data(void);
    void               release_data(void);
    void               alloc_nulval(const unsigned char* value);
    void               init_data(void);
    void*              ptr_data(void) { return m_data; }
    void*              ptr_nulval(void) { return m_nulval; }

    // Private data area
    unsigned char* m_data;       //!< Data vector
    unsigned char* m_nulval;     //!< NULL value
};

#endif /* GFITSTABLEBYTECOL_HPP */
