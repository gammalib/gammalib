/***************************************************************************
 *          GFitsTableLongCol.hpp  - FITS table long column class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsTableLongCol.hpp
 * @brief GFitsTableLongCol class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSTABLELONGCOL_HPP
#define GFITSTABLELONGCOL_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsTableCol.hpp"


/***********************************************************************//**
 * @class GFitsTableLongCol
 *
 * @brief Interface for FITS table long integer column
 *
 * This class implements a FITS table long integer column.
 ***************************************************************************/
class GFitsTableLongCol : public GFitsTableCol {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GFitsTableLongCol& column);

public:
    // Constructors and destructors
    GFitsTableLongCol(void);
    GFitsTableLongCol(const std::string& name, const int& length,
                      const int& size = 1);
    GFitsTableLongCol(const GFitsTableLongCol& column);
    virtual ~GFitsTableLongCol(void);

    // Operators
    GFitsTableLongCol& operator= (const GFitsTableLongCol& column);
    long&              operator() (const int& row, const int& inx = 0);
    const long&        operator() (const int& row, const int& inx = 0) const;

    // Methods
    std::string string(const int& row, const int& col = 0);
    double      real(const int& row, const int& col = 0);
    int         integer(const int& row, const int& col = 0);
    long*       data(void) { return m_data; }
    void        nullval(const long* value);
    long*       nullval(void) { return m_nulval; }

private:
    // Private methods
    void               init_members(void);
    void               copy_members(const GFitsTableLongCol& column);
    void               free_members(void);
    GFitsTableLongCol* clone(void) const;
    std::string        ascii_format(void) const;
    std::string        binary_format(void) const;
    void               alloc_data(void);
    void               release_data(void);
    void               alloc_nulval(const long* value);
    void               init_data(void);
    void*              ptr_data(void) { return m_data; }
    void*              ptr_nulval(void) { return m_nulval; }

    // Private data area
    long* m_data;       //!< Data vector
    long* m_nulval;     //!< NULL value
};

#endif /* GFITSTABLELONGCOL_HPP */
