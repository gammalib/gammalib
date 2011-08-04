/***************************************************************************
 *    GFitsTableUShortCol.hpp  - FITS table unsigned short column class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @file GFitsTableUShortCol.hpp
 * @brief FITS table unsigned short integer column class interface definition
 * @author J. Knodlseder
 */

#ifndef GFITSTABLEUSHORTCOL_HPP
#define GFITSTABLEUSHORTCOL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GFitsTableCol.hpp"


/***********************************************************************//**
 * @class GFitsTableUShortCol
 *
 * @brief FITS table unsigned short integer column
 *
 * This class implements a FITS table unsigned short integer column.
 ***************************************************************************/
class GFitsTableUShortCol : public GFitsTableCol {

public:
    // Constructors and destructors
    GFitsTableUShortCol(void);
    GFitsTableUShortCol(const std::string& name, const int& length,
                        const int& size = 1);
    GFitsTableUShortCol(const GFitsTableUShortCol& column);
    virtual ~GFitsTableUShortCol(void);

    // Operators
    GFitsTableUShortCol&  operator= (const GFitsTableUShortCol& column);
    unsigned short&       operator() (const int& row, const int& inx = 0);
    const unsigned short& operator() (const int& row, const int& inx = 0) const;

    // Implement virtual methods
    virtual std::string string(const int& row, const int& col = 0);
    virtual double      real(const int& row, const int& col = 0);
    virtual int         integer(const int& row, const int& col = 0);
    virtual void        insert(const int& rownum, const int& nrows);
    virtual void        remove(const int& rownum, const int& nrows);
    
    // Other methods
    unsigned short* data(void) { return m_data; }
    void            nulval(const unsigned short* value);
    unsigned short* nulval(void) { return m_nulval; }

private:
    // Private methods
    void                 init_members(void);
    void                 copy_members(const GFitsTableUShortCol& column);
    void                 free_members(void);
    GFitsTableUShortCol* clone(void) const;
    std::string          ascii_format(void) const;
    std::string          binary_format(void) const;
    void                 alloc_data(void);
    void                 release_data(void);
    void                 alloc_nulval(const unsigned short* value);
    void                 init_data(void);
    void*                ptr_data(void) { return m_data; }
    void*                ptr_nulval(void) { return m_nulval; }

    // Private data area
    unsigned short* m_data;       //!< Data vector
    unsigned short* m_nulval;     //!< NULL value
};

#endif /* GFITSTABLEUSHORTCOL_HPP */
