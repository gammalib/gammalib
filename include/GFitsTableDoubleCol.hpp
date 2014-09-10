/***************************************************************************
 *         GFitsTableDoubleCol.hpp - FITS table double column class        *
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
 * @file GFitsTableDoubleCol.hpp
 * @brief FITS table double column class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GFITSTABLEDOUBLECOL_HPP
#define GFITSTABLEDOUBLECOL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GFitsTableCol.hpp"


/***********************************************************************//**
 * @class GFitsTableDoubleCol
 *
 * @brief FITS table double column
 *
 * This class implements a FITS table double column.
 ***************************************************************************/
class GFitsTableDoubleCol : public GFitsTableCol {

public:
    // Constructors and destructors
    GFitsTableDoubleCol(void);
    GFitsTableDoubleCol(const std::string& name, const int& length,
                        const int& size = 1);
    GFitsTableDoubleCol(const GFitsTableDoubleCol& column);
    virtual ~GFitsTableDoubleCol(void);

    // Operators
    GFitsTableDoubleCol& operator=(const GFitsTableDoubleCol& column);
    double&              operator()(const int& row, const int& inx = 0);
    const double&        operator()(const int& row, const int& inx = 0) const;

    // Implement virtual methods
    virtual void                 clear(void);
    virtual GFitsTableDoubleCol* clone(void) const;
    virtual std::string          classname(void) const;
    virtual std::string          string(const int& row, const int& col = 0) const;
    virtual double               real(const int& row, const int& col = 0) const;
    virtual int                  integer(const int& row, const int& col = 0) const;
    virtual void                 insert(const int& row, const int& nrows);
    virtual void                 remove(const int& row, const int& nrows);
    virtual bool                 is_loaded(void) const;
    
    // Other methods
    double* data(void);
    double* nulval(void);
    void    nulval(const double* value);

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsTableDoubleCol& column);
    void free_members(void);
    void alloc_nulval(const double* value);

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
    double* m_data;       //!< Data vector
    double* m_nulval;     //!< NULL value
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GFitsTableDoubleCol").
 ***************************************************************************/
inline
std::string GFitsTableDoubleCol::classname(void) const
{
    return ("GFitsTableDoubleCol");
}


/***********************************************************************//**
 * @brief Checks if column has been loaded
 *
 * @return True if column has been loaded, false otherwise
 ***************************************************************************/
inline
bool GFitsTableDoubleCol::is_loaded(void) const
{
    return (m_data != NULL);
}


/***********************************************************************//**
 * @brief Returns pointer to column data
 *
 * @return Pointer to column data.
 ***************************************************************************/
inline
double* GFitsTableDoubleCol::data(void)
{
    return m_data;
}


/***********************************************************************//**
 * @brief Returns pointer to nul value
 *
 * @return Pointer to nul value.
 ***************************************************************************/
inline
double* GFitsTableDoubleCol::nulval(void)
{
    return m_nulval;
}


/***********************************************************************//**
 * @brief Returns void pointer to column data
 *
 * @return Void pointer to column data.
 ***************************************************************************/
inline
void* GFitsTableDoubleCol::ptr_data(const int& index)
{
    return (m_data+index);
}


/***********************************************************************//**
 * @brief Returns void pointer to nul value
 *
 * @return Void pointer to nul value.
 ***************************************************************************/
inline
void* GFitsTableDoubleCol::ptr_nulval(void)
{
    return m_nulval;
}

#endif /* GFITSTABLEDOUBLECOL_HPP */
