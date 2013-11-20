/***************************************************************************
 *   GFitsTableCFloatCol.hpp - FITS table single precision complex column  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2013 by Juergen Knoedlseder                         *
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
 * @file GFitsTableCFloatCol.hpp
 * @brief FITS table float complex column class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GFITSTABLECFLOATCOL_HPP
#define GFITSTABLECFLOATCOL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GFits.hpp"
#include "GFitsTableCol.hpp"


/***********************************************************************//**
 * @class GFitsTableCFloatCol
 *
 * @brief FITS table float complex column
 *
 * This class implements a FITS table float complex column.
 ***************************************************************************/
class GFitsTableCFloatCol : public GFitsTableCol {

public:
    // Constructors and destructors
    GFitsTableCFloatCol(void);
    GFitsTableCFloatCol(const std::string& name, const int& length,
                         const int& size = 1);
    GFitsTableCFloatCol(const GFitsTableCFloatCol& column);
    virtual ~GFitsTableCFloatCol(void);

    // Operators
    GFitsTableCFloatCol& operator=(const GFitsTableCFloatCol& column);
    GFits::cfloat&        operator()(const int& row, const int& inx = 0);
    const GFits::cfloat&  operator()(const int& row, const int& inx = 0) const;

    // Implement virtual methods
    virtual void                 clear(void);
    virtual GFitsTableCFloatCol* clone(void) const;
    virtual std::string          string(const int& row, const int& col = 0) const;
    virtual double               real(const int& row, const int& col = 0) const;
    virtual int                  integer(const int& row, const int& col = 0) const;
    virtual void                 insert(const int& row, const int& nrows);
    virtual void                 remove(const int& row, const int& nrows);
    
    // Other methods
    GFits::cfloat* data(void) { return m_data; }
    void           nulval(const GFits::cfloat* value);
    GFits::cfloat* nulval(void) { return m_nulval; }

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsTableCFloatCol& column);
    void free_members(void);
    void alloc_nulval(const GFits::cfloat* value);

    // Implemented virtual base class methods
    virtual void        alloc_data(void);
    virtual void        init_data(void);
    virtual void        fetch_data(void) const;
    virtual void        copy_data(const GFitsTableCol& column);
    virtual void        release_data(void);
    virtual void*       ptr_data(const int& index = 0) { return m_data+index; }
    virtual void*       ptr_nulval(void) { return m_nulval; }
    virtual std::string ascii_format(void) const;

    // Private data area
    GFits::cfloat* m_data;       //!< Data vector
    GFits::cfloat* m_nulval;     //!< NULL value
};

#endif /* GFITSTABLECFLOATCOL_HPP */
