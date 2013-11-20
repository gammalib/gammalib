/***************************************************************************
 *         GFitsTableByteCol.hpp  - FITS table Byte column class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2012 by Juergen Knoedlseder                         *
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
 * @file GFitsTableByteCol.hpp
 * @brief FITS table Byte column class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GFITSTABLEBYTECOL_HPP
#define GFITSTABLEBYTECOL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GFitsTableCol.hpp"


/***********************************************************************//**
 * @class GFitsTableByteCol
 *
 * @brief FITS table Byte column
 *
 * This class implements a FITS table Byte column.
 ***************************************************************************/
class GFitsTableByteCol : public GFitsTableCol {

public:
    // Constructors and destructors
    GFitsTableByteCol(void);
    GFitsTableByteCol(const std::string& name, const int& length,
                      const int& size = 1);
    GFitsTableByteCol(const GFitsTableByteCol& column);
    virtual ~GFitsTableByteCol(void);

    // Operators
    GFitsTableByteCol&   operator= (const GFitsTableByteCol& column);
    unsigned char&       operator() (const int& row, const int& inx = 0);
    const unsigned char& operator() (const int& row, const int& inx = 0) const;

    // Implement virtual methods
    virtual void               clear(void);
    virtual GFitsTableByteCol* clone(void) const;
    virtual std::string        string(const int& row, const int& col = 0) const;
    virtual double             real(const int& row, const int& col = 0) const;
    virtual int                integer(const int& row, const int& col = 0) const;
    virtual void               insert(const int& rownum, const int& nrows);
    virtual void               remove(const int& rownum, const int& nrows);
    
    // Other methods
    unsigned char* data(void) { return m_data; }
    void           nulval(const unsigned char* value);
    unsigned char* nulval(void) { return m_nulval; }

private:
    // Private methods
    void        init_members(void);
    void        copy_members(const GFitsTableByteCol& column);
    void        free_members(void);
    std::string ascii_format(void) const;
    std::string binary_format(void) const;
    void        alloc_data(void);
    void        release_data(void);
    void        alloc_nulval(const unsigned char* value);
    void        init_data(void);
    void*       ptr_data(const int& index = 0) { return m_data+index; }
    void*       ptr_nulval(void) { return m_nulval; }

    // Private data area
    unsigned char* m_data;       //!< Data vector
    unsigned char* m_nulval;     //!< NULL value
};

#endif /* GFITSTABLEBYTECOL_HPP */
