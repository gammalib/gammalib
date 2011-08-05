/***************************************************************************
 *         GFitsTableFloatCol.hpp  - FITS table float column class         *
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
 * @file GFitsTableFloatCol.hpp
 * @brief FITS table float column class interface definition
 * @author J. Knodlseder
 */

#ifndef GFITSTABLEFLOATCOL_HPP
#define GFITSTABLEFLOATCOL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GFitsTableCol.hpp"


/***********************************************************************//**
 * @class GFitsTableFloatCol
 *
 * @brief FITS table float column
 *
 * This class implements a FITS table float column.
 ***************************************************************************/
class GFitsTableFloatCol : public GFitsTableCol {

public:
    // Constructors and destructors
    GFitsTableFloatCol(void);
    GFitsTableFloatCol(const std::string& name, const int& length,
                       const int& size = 1);
    GFitsTableFloatCol(const GFitsTableFloatCol& column);
    virtual ~GFitsTableFloatCol(void);

    // Operators
    GFitsTableFloatCol& operator= (const GFitsTableFloatCol& column);
    float&              operator() (const int& row, const int& inx = 0);
    const float&        operator() (const int& row, const int& inx = 0) const;

    // Implement virtual methods
    virtual std::string string(const int& row, const int& col = 0) const;
    virtual double      real(const int& row, const int& col = 0) const;
    virtual int         integer(const int& row, const int& col = 0) const;
    virtual void        insert(const int& rownum, const int& nrows);
    virtual void        remove(const int& rownum, const int& nrows);
    
    // Other methods
    float* data(void) { return m_data; }
    void   nulval(const float* value);
    float* nulval(void) { return m_nulval; }

private:
    // Private methods
    void                init_members(void);
    void                copy_members(const GFitsTableFloatCol& column);
    void                free_members(void);
    GFitsTableFloatCol* clone(void) const;
    std::string         ascii_format(void) const;
    std::string         binary_format(void) const;
    void                alloc_data(void);
    void                release_data(void);
    void                alloc_nulval(const float* value);
    void                init_data(void);
    void*               ptr_data(void) { return m_data; }
    void*               ptr_nulval(void) { return m_nulval; }

    // Private data area
    float* m_data;       //!< Data vector
    float* m_nulval;     //!< NULL value
};

#endif /* GFITSTABLEFLOATCOL_HPP */
