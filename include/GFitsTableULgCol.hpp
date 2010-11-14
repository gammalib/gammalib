/***************************************************************************
 *      GFitsTableULgCol.hpp  - FITS table unsigned long column class      *
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
 * @file GFitsTableULgCol.hpp
 * @brief GFitsTableULgCol class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSTABLEULGCOL_HPP
#define GFITSTABLEULGCOL_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsTableCol.hpp"


/***********************************************************************//**
 * @class GFitsTableULgCol
 *
 * @brief Interface for FITS table unsigned long integer column
 *
 * This class implements a FITS table unsigned long integer column.
 ***************************************************************************/
class GFitsTableULgCol : public GFitsTableCol {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GFitsTableULgCol& column);

public:
    // Constructors and destructors
    GFitsTableULgCol();
    GFitsTableULgCol(const std::string& name, const int& length,
                     const int& size = 1);
    GFitsTableULgCol(const GFitsTableULgCol& column);
    virtual ~GFitsTableULgCol();

    // Operators
    GFitsTableULgCol&    operator= (const GFitsTableULgCol& column);
    unsigned long&       operator() (const int& row, const int& inx = 0);
    const unsigned long& operator() (const int& row, const int& inx = 0) const;

    // Methods
    std::string    string(const int& row, const int& col = 0);
    double         real(const int& row, const int& col = 0);
    int            integer(const int& row, const int& col = 0);
    unsigned long* data(void);
    void           set_nullval(const unsigned long* value);

private:
    // Private methods
    void              init_members(void);
    void              copy_members(const GFitsTableULgCol& column);
    void              free_members(void);
    void              save(void);
    GFitsTableULgCol* clone(void) const;
    std::string       ascii_format(void) const;
    std::string       binary_format(void) const;
    void              alloc_data(void);
    void              init_data(void);
    void              fetch_data(void);
    void*             ptr_data(void) { return m_data; }
    void*             ptr_nulval(void) { return m_nulval; }

    // Private data area
    unsigned long* m_data;          //!< Data area
    unsigned long* m_nulval;        //!< NULL value
};

#endif /* GFITSTABLEULGCOL_HPP */
