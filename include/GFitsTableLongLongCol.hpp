/***************************************************************************
 *      GFitsTableLongLongCol.hpp - FITS table long long column class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
 * @file GFitsTableLongLongCol.hpp
 * @brief FITS table long long integer column class interface definition
 * @author Juergen Knodlseder
 */

#ifndef GFITSTABLELONGLONGCOL_HPP
#define GFITSTABLELONGLONGCOL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GFitsTableCol.hpp"


/***********************************************************************//**
 * @class GFitsTableLongLongCol
 *
 * @brief FITS table long long integer column
 *
 * This class implements a FITS table long long integer column.
 ***************************************************************************/
class GFitsTableLongLongCol : public GFitsTableCol {

public:
    // Constructors and destructors
    GFitsTableLongLongCol(void);
    GFitsTableLongLongCol(const std::string& name, const int& length,
                          const int& size = 1);
    GFitsTableLongLongCol(const GFitsTableLongLongCol& column);
    virtual ~GFitsTableLongLongCol(void);

    // Operators
    GFitsTableLongLongCol& operator=(const GFitsTableLongLongCol& column);
    long long&             operator()(const int& row, const int& inx = 0);
    const long long&       operator()(const int& row, const int& inx = 0) const;

    // Implement virtual methods
    virtual void                   clear(void);
    virtual GFitsTableLongLongCol* clone(void) const;
    virtual std::string            classname(void) const;
    virtual std::string            string(const int& row, const int& col = 0) const;
    virtual double                 real(const int& row, const int& col = 0) const;
    virtual int                    integer(const int& row, const int& col = 0) const;
    virtual void                   insert(const int& row, const int& nrows);
    virtual void                   remove(const int& row, const int& nrows);
    virtual bool                   is_loaded(void) const;
    
    // Other methods
    long long*  data(void);
    long long*  nulval(void);
    void        nulval(const long long* value);

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsTableLongLongCol& column);
    void free_members(void);
    void alloc_nulval(const long long* value);

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
    long long* m_data;       //!< Data vector
    long long* m_nulval;     //!< NULL value
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GFitsTableLongLongCol").
 ***************************************************************************/
inline
std::string GFitsTableLongLongCol::classname(void) const
{
    return ("GFitsTableLongLongCol");
}


/***********************************************************************//**
 * @brief Checks if column has been loaded
 *
 * @return True if column has been loaded, false otherwise
 ***************************************************************************/
inline
bool GFitsTableLongLongCol::is_loaded(void) const
{
    return (m_data != NULL);
}


/***********************************************************************//**
 * @brief Returns pointer to column data
 *
 * @return Pointer to column data.
 ***************************************************************************/
inline
long long* GFitsTableLongLongCol::data(void)
{
    return m_data;
}


/***********************************************************************//**
 * @brief Returns pointer to nul value
 *
 * @return Pointer to nul value.
 ***************************************************************************/
inline
long long* GFitsTableLongLongCol::nulval(void)
{
    return m_nulval;
}


/***********************************************************************//**
 * @brief Returns void pointer to column data
 *
 * @return Void pointer to column data.
 ***************************************************************************/
inline
void* GFitsTableLongLongCol::ptr_data(const int& index)
{
    return (m_data+index);
}


/***********************************************************************//**
 * @brief Returns void pointer to nul value
 *
 * @return Void pointer to nul value.
 ***************************************************************************/
inline
void* GFitsTableLongLongCol::ptr_nulval(void)
{
    return m_nulval;
}

#endif /* GFITSTABLELONGLONGCOL_HPP */
