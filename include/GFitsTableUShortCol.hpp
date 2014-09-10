/***************************************************************************
 *     GFitsTableUShortCol.hpp - FITS table unsigned short column class    *
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
 * @file GFitsTableUShortCol.hpp
 * @brief FITS table unsigned short integer column class interface definition
 * @author Juergen Knoedlseder
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
    GFitsTableUShortCol&  operator=(const GFitsTableUShortCol& column);
    unsigned short&       operator()(const int& row, const int& inx = 0);
    const unsigned short& operator()(const int& row, const int& inx = 0) const;

    // Implement virtual methods
    virtual void                 clear(void);
    virtual GFitsTableUShortCol* clone(void) const;
    virtual std::string          classname(void) const;
    virtual std::string          string(const int& row, const int& col = 0) const;
    virtual double               real(const int& row, const int& col = 0) const;
    virtual int                  integer(const int& row, const int& col = 0) const;
    virtual void                 insert(const int& row, const int& nrows);
    virtual void                 remove(const int& row, const int& nrows);
    virtual bool                 is_loaded(void) const;
    
    // Other methods
    unsigned short* data(void);
    unsigned short* nulval(void);
    void            nulval(const unsigned short* value);

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsTableUShortCol& column);
    void free_members(void);
    void alloc_nulval(const unsigned short* value);

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
    unsigned short* m_data;       //!< Data vector
    unsigned short* m_nulval;     //!< NULL value
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GFitsTableUShortCol").
 ***************************************************************************/
inline
std::string GFitsTableUShortCol::classname(void) const
{
    return ("GFitsTableUShortCol");
}


/***********************************************************************//**
 * @brief Checks if column has been loaded
 *
 * @return True if column has been loaded, false otherwise
 ***************************************************************************/
inline
bool GFitsTableUShortCol::is_loaded(void) const
{
    return (m_data != NULL);
}


/***********************************************************************//**
 * @brief Returns pointer to column data
 *
 * @return Pointer to column data.
 ***************************************************************************/
inline
unsigned short* GFitsTableUShortCol::data(void)
{
    return m_data;
}


/***********************************************************************//**
 * @brief Returns pointer to nul value
 *
 * @return Pointer to nul value.
 ***************************************************************************/
inline
unsigned short* GFitsTableUShortCol::nulval(void)
{
    return m_nulval;
}


/***********************************************************************//**
 * @brief Returns void pointer to column data
 *
 * @return Void pointer to column data.
 ***************************************************************************/
inline
void* GFitsTableUShortCol::ptr_data(const int& index)
{
    return (m_data+index);
}


/***********************************************************************//**
 * @brief Returns void pointer to nul value
 *
 * @return Void pointer to nul value.
 ***************************************************************************/
inline
void* GFitsTableUShortCol::ptr_nulval(void)
{
    return m_nulval;
}

#endif /* GFITSTABLEUSHORTCOL_HPP */
