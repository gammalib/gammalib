/***************************************************************************
 *        GFitsTableBoolCol.hpp  - FITS table boolean column class         *
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
 * @file GFitsTableBoolCol.hpp
 * @brief FITS table Boolean column class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GFITSTABLEBOOLCOL_HPP
#define GFITSTABLEBOOLCOL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GFitsTableCol.hpp"


/***********************************************************************//**
 * @class GFitsTableBoolCol
 *
 * @brief FITS table Boolean column
 *
 * This class implements a FITS table Boolean column.
 *
 * @todo Each Bit is actually stored in one Byte. To save memory a more
 *       compact storage scheme should be implemented.
 ***************************************************************************/
class GFitsTableBoolCol : public GFitsTableCol {

public:
    // Constructors and destructors
    GFitsTableBoolCol(void);
    GFitsTableBoolCol(const std::string& name, const int& length,
                      const int& size = 1);
    GFitsTableBoolCol(const GFitsTableBoolCol& column);
    virtual ~GFitsTableBoolCol(void);

    // Operators
    GFitsTableBoolCol& operator= (const GFitsTableBoolCol& column);
    bool&              operator() (const int& row, const int& inx = 0);
    const bool&        operator() (const int& row, const int& inx = 0) const;

    // Implement virtual methods
    virtual void               clear(void);
    virtual GFitsTableBoolCol* clone(void) const;
    virtual std::string        string(const int& row, const int& col = 0) const;
    virtual double             real(const int& row, const int& col = 0) const;
    virtual int                integer(const int& row, const int& col = 0) const;
    virtual void               insert(const int& rownum, const int& nrows);
    virtual void               remove(const int& rownum, const int& nrows);
    
    // Other methods
    bool* data(void) { return m_data; }
    void  nulval(const bool* value);
    bool* nulval(void) { return m_nulval; }

private:
    // Private methods
    void        init_members(void);
    void        copy_members(const GFitsTableBoolCol& column);
    void        free_members(void);
    void        save(void);
    std::string ascii_format(void) const;
    std::string binary_format(void) const;
    void        alloc_data(void);
    void        release_data(void);
    void        alloc_nulval(const bool* value);
    void        init_data(void);
    void        fetch_data(void) const;
    void        alloc_buffer(void);
    void        free_buffer(void);
    void*       ptr_data(void) { return m_buffer; }
    void*       ptr_nulval(void) { return m_nulval; }

    // Private data area
    bool* m_data;       //!< Data area
    char* m_buffer;     //!< Data area for CFITSIO transfer
    bool* m_nulval;     //!< NULL value
};

#endif /* GFITSTABLEBOOLCOL_HPP */
