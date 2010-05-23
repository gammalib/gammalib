/***************************************************************************
 *        GFitsTableLlgCol.hpp  - FITS table long long column class        *
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
 * @file GFitsTableLlgCol.hpp
 * @brief GFitsTableLlgCol class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSTABLELLGCOL_HPP
#define GFITSTABLELLGCOL_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsTableCol.hpp"


/***********************************************************************//**
 * @class GFitsTableLlgCol
 *
 * @brief Interface for FITS table long long integer column
 *
 * This class implements a FITS table long long integer column.
 *
 * TODO: integer() method returns actually a int which may lead to an
 * overflow. Better return long long.
 ***************************************************************************/
class GFitsTableLlgCol : public GFitsTableCol {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GFitsTableLlgCol& column);

public:
    // Constructors and destructors
    GFitsTableLlgCol(void);
    GFitsTableLlgCol(const std::string& name, const int& length,
                     const int& size = 1);
    GFitsTableLlgCol(const GFitsTableLlgCol& column);
    virtual ~GFitsTableLlgCol(void);

    // Operators
    GFitsTableLlgCol& operator= (const GFitsTableLlgCol& column);
    long long&        operator() (const int& row, const int& inx = 0);
    const long long&  operator() (const int& row, const int& inx = 0) const;

    // Methods
    std::string string(const int& row, const int& col = 0);
    double      real(const int& row, const int& col = 0);
    int         integer(const int& row, const int& col = 0);
    long long*  data(void);
    void        set_nullval(const long long* value);

private:
    // Private methods
    void              init_members(void);
    void              copy_members(const GFitsTableLlgCol& column);
    void              free_members(void);
    void              save(void);
    GFitsTableLlgCol* clone(void) const;
    std::string       ascii_format(void) const;
    std::string       binary_format(void) const;
    void              alloc_data(void);
    void              init_data(void);
    void              fetch_data(void);
    void*             ptr_data(void) { return m_data; }
    void*             ptr_nulval(void) { return m_nulval; }

    // Private data area
    long long* m_data;          //!< Data area
    long long* m_nulval;        //!< NULL value
};

#endif /* GFITSTABLELLGCOL_HPP */
