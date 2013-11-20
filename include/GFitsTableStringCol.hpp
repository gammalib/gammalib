/***************************************************************************
 *        GFitsTableStringCol.hpp  - FITS table string column class        *
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
 * @file GFitsTableStringCol.hpp
 * @brief FITS table string column class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GFITSTABLESTRINGCOL_HPP
#define GFITSTABLESTRINGCOL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GFitsTableCol.hpp"


/***********************************************************************//**
 * @class GFitsTableStringCol
 *
 * @brief FITS table string column
 *
 * This class implements a FITS table string column.
 ***************************************************************************/
class GFitsTableStringCol : public GFitsTableCol {

public:
    // Constructors and destructors
    GFitsTableStringCol(void);
    GFitsTableStringCol(const std::string& name, const int& length,
                        const int& width, const int& size = 1);
    GFitsTableStringCol(const GFitsTableStringCol& column);
    virtual ~GFitsTableStringCol(void);

    // Operators
    GFitsTableStringCol&  operator= (const GFitsTableStringCol& column);
    std::string&          operator() (const int& row, const int& inx = 0);
    const std::string&    operator() (const int& row, const int& inx = 0) const;

    // Implement virtual methods
    virtual void                 clear(void);
    virtual GFitsTableStringCol* clone(void) const;
    virtual std::string          string(const int& row, const int& col = 0) const;
    virtual double               real(const int& row, const int& col = 0) const;
    virtual int                  integer(const int& row, const int& col = 0) const;
    virtual void                 insert(const int& rownum, const int& nrows);
    virtual void                 remove(const int& rownum, const int& nrows);
    
    // Other methods
    std::string* data(void) { return m_data; }
    void         nulval(const std::string& value);
    char*        nulval(void) { return m_nulval; }

private:
    // Private methods
    void                 init_members(void);
    void                 copy_members(const GFitsTableStringCol& column);
    void                 free_members(void);
    void                 save(void);
    std::string          ascii_format(void) const;
    std::string          binary_format(void) const;
    void                 alloc_data(void);
    void                 release_data(void);
    void                 alloc_nulval(const std::string& value);
    void                 init_data(void);
    void                 fetch_data(void) const;
    void                 alloc_buffer(void) const;
    void                 free_buffer(void) const;
    void*                ptr_data(const int& index = 0) { return m_buffer+index; }
    void*                ptr_nulval(void) { return m_nulval; }

    // Private data area
    std::string*    m_data;    //!< Data area
    mutable char**  m_buffer;  //!< Data area for CFITSIO transfer
    char*           m_nulval;  //!< NULL string
};

#endif /* GFITSTABLESTRINGCOL_HPP */
