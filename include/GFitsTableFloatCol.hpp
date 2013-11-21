/***************************************************************************
 *          GFitsTableFloatCol.hpp - FITS table float column class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
    GFitsTableFloatCol& operator=(const GFitsTableFloatCol& column);
    float&              operator()(const int& row, const int& inx = 0);
    const float&        operator()(const int& row, const int& inx = 0) const;

    // Implement virtual methods
    virtual void                clear(void);
    virtual GFitsTableFloatCol* clone(void) const;
    virtual std::string         string(const int& row, const int& col = 0) const;
    virtual double              real(const int& row, const int& col = 0) const;
    virtual int                 integer(const int& row, const int& col = 0) const;
    virtual void                insert(const int& row, const int& nrows);
    virtual void                remove(const int& row, const int& nrows);
    virtual bool                isloaded(void) const;
    
    // Other methods
    float* data(void);
    float* nulval(void);
    void   nulval(const float* value);

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsTableFloatCol& column);
    void free_members(void);
    void alloc_nulval(const float* value);

    // Implemented virtual base class methods
    virtual void        alloc_data(void);
    virtual void        init_data(void);
    virtual void        fetch_data(void) const;
    virtual void        release_data(void);
    virtual void*       ptr_data(const int& index = 0);
    virtual void*       ptr_nulval(void);
    virtual std::string ascii_format(void) const;

    // Private data area
    float* m_data;       //!< Data vector
    float* m_nulval;     //!< NULL value
};


/***********************************************************************//**
 * @brief Checks if column has been loaded
 *
 * @return True if column has been loaded, false otherwise
 ***************************************************************************/
inline
bool GFitsTableFloatCol::isloaded(void) const
{
    return (m_data != NULL);
}


/***********************************************************************//**
 * @brief Returns pointer to column data
 *
 * @return Pointer to column data.
 ***************************************************************************/
inline
float* GFitsTableFloatCol::data(void)
{
    return m_data;
}


/***********************************************************************//**
 * @brief Returns pointer to nul value
 *
 * @return Pointer to nul value.
 ***************************************************************************/
inline
float* GFitsTableFloatCol::nulval(void)
{
    return m_nulval;
}


/***********************************************************************//**
 * @brief Returns void pointer to column data
 *
 * @return Void pointer to column data.
 ***************************************************************************/
inline
void* GFitsTableFloatCol::ptr_data(const int& index)
{
    return (m_data+index);
}


/***********************************************************************//**
 * @brief Returns void pointer to nul value
 *
 * @return Void pointer to nul value.
 ***************************************************************************/
inline
void* GFitsTableFloatCol::ptr_nulval(void)
{
    return m_nulval;
}

#endif /* GFITSTABLEFLOATCOL_HPP */
