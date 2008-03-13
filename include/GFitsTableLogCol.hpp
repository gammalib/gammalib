/***************************************************************************
 *         GFitsTableLogCol.hpp  - FITS table logical column class         *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsTableLogCol.hpp
 * @brief GFitsTableLogCol class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSTABLELOGCOL_HPP
#define GFITSTABLELOGCOL_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsCfitsio.hpp"
#include "GFitsTableCol.hpp"

/* __ Namespaces _________________________________________________________ */


/* __ Structures _________________________________________________________ */


/***********************************************************************//**
 * @class GFitsTableLogCol
 *
 * @brief Interface for FITS table logical column
 *
 * This class implements a FITS table logical column.
 ***************************************************************************/
class GFitsTableLogCol : public GFitsTableCol {

    // I/O friends
    friend ostream& operator<< (ostream& os, const GFitsTableLogCol& column);

public:
    // Constructors and destructors
    GFitsTableLogCol();
    GFitsTableLogCol(const std::string& name, const int& length,
                     const int& size = 1);
    GFitsTableLogCol(const GFitsTableLogCol& column);
    virtual ~GFitsTableLogCol();

    // Operators
    GFitsTableLogCol& operator= (const GFitsTableLogCol& column);
    char&             operator() (const int& row, const int& inx = 0);
    const char&       operator() (const int& row, const int& inx = 0) const;

    // Methods
    std::string string(const int& row, const int& col = 0);
    double      real(const int& row, const int& col = 0);
    int         integer(const int& row, const int& col = 0);
    char*       data(void);
    void        set_nullval(const char* value);

private:
    // Private methods
    void              init_members(void);
    void              copy_members(const GFitsTableLogCol& column);
    void              free_members(void);
    void              save(void);
    GFitsTableLogCol* clone(void) const;
    std::string       ascii_format(void) const;
    std::string       binary_format(void) const;
    void              alloc_data(void);
    void              init_data(void);
    void              fetch_data(void);
    void*             ptr_data(void) { return m_data; }
    void*             ptr_nulval(void) { return m_nulval; }

    // Private data area
    char* m_data;       //!< Data area
    char* m_nulval;     //!< NULL value
};

#endif /* GFITSTABLELOGCOL_HPP */
