/***************************************************************************
 *          GFitsTableLongCol.hpp - FITS table long column class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Jurgen Knodlseder                           *
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
 * @file GFitsTableLongCol.hpp
 * @brief FITS table long integer column class interface definition
 * @author J. Knodlseder
 */

#ifndef GFITSTABLELONGCOL_HPP
#define GFITSTABLELONGCOL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GFitsTableCol.hpp"


/***********************************************************************//**
 * @class GFitsTableLongCol
 *
 * @brief FITS table long integer column
 *
 * This class implements a FITS table long integer column.
 ***************************************************************************/
class GFitsTableLongCol : public GFitsTableCol {

public:
    // Constructors and destructors
    GFitsTableLongCol(void);
    GFitsTableLongCol(const std::string& name, const int& length,
                      const int& size = 1);
    GFitsTableLongCol(const GFitsTableLongCol& column);
    virtual ~GFitsTableLongCol(void);

    // Operators
    GFitsTableLongCol& operator=(const GFitsTableLongCol& column);
    long&              operator()(const int& row, const int& inx = 0);
    const long&        operator()(const int& row, const int& inx = 0) const;

    // Implement virtual methods
    virtual void               clear(void);
    virtual GFitsTableLongCol* clone(void) const;
    virtual std::string        classname(void) const;
    virtual std::string        string(const int& row, const int& col = 0) const;
    virtual double             real(const int& row, const int& col = 0) const;
    virtual int                integer(const int& row, const int& col = 0) const;
    virtual void               insert(const int& row, const int& nrows);
    virtual void               remove(const int& row, const int& nrows);
    virtual bool               is_loaded(void) const;
    
    // Other methods
    long* data(void);
    long* nulval(void);
    void  nulval(const long* value);

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsTableLongCol& column);
    void free_members(void);
    void alloc_nulval(const long* value);

    // Implemented virtual base class methods
    virtual void        alloc_data(void);
    virtual void        init_data(void);
    virtual void        fetch_data(void) const;
    virtual void        resize_data(const int& index, const int& number);
    virtual void        release_data(void);
    virtual void*       ptr_data(const int& index = 0);
    virtual void*       ptr_nulval(void);
    virtual std::string ascii_format(void) const;

    // Private data area
    long* m_data;       //!< Data vector
    long* m_nulval;     //!< NULL value
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GFitsTableLongCol").
 ***************************************************************************/
inline
std::string GFitsTableLongCol::classname(void) const
{
    return ("GFitsTableLongCol");
}


/***********************************************************************//**
 * @brief Checks if column has been loaded
 *
 * @return True if column has been loaded, false otherwise
 ***************************************************************************/
inline
bool GFitsTableLongCol::is_loaded(void) const
{
    return (m_data != NULL);
}


/***********************************************************************//**
 * @brief Returns pointer to column data
 *
 * @return Pointer to column data.
 ***************************************************************************/
inline
long* GFitsTableLongCol::data(void)
{
    return m_data;
}


/***********************************************************************//**
 * @brief Returns pointer to nul value
 *
 * @return Pointer to nul value.
 ***************************************************************************/
inline
long* GFitsTableLongCol::nulval(void)
{
    return m_nulval;
}


/***********************************************************************//**
 * @brief Returns void pointer to column data
 *
 * @return Void pointer to column data.
 ***************************************************************************/
inline
void* GFitsTableLongCol::ptr_data(const int& index)
{
    return (m_data+index);
}


/***********************************************************************//**
 * @brief Returns void pointer to nul value
 *
 * @return Void pointer to nul value.
 ***************************************************************************/
inline
void* GFitsTableLongCol::ptr_nulval(void)
{
    return m_nulval;
}

#endif /* GFITSTABLELONGCOL_HPP */
