/***************************************************************************
 *         GFitsTableShortCol.hpp  - FITS table short column class         *
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
 * @file GFitsTableShortCol.hpp
 * @brief GFitsTableShortCol class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSTABLESHORTCOL_HPP
#define GFITSTABLESHORTCOL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GFitsTableCol.hpp"


/***********************************************************************//**
 * @class GFitsTableShortCol
 *
 * @brief Interface for FITS table short integer column
 *
 * This class implements a FITS table short integer column.
 ***************************************************************************/
class GFitsTableShortCol : public GFitsTableCol {

public:
    // Constructors and destructors
    GFitsTableShortCol(void);
    GFitsTableShortCol(const std::string& name, const int& length,
                       const int& size = 1);
    GFitsTableShortCol(const GFitsTableShortCol& column);
    virtual ~GFitsTableShortCol(void);

    // Operators
    GFitsTableShortCol& operator= (const GFitsTableShortCol& column);
    short&              operator() (const int& row, const int& inx = 0);
    const short&        operator() (const int& row, const int& inx = 0) const;

    // Methods
    std::string string(const int& row, const int& col = 0);
    double      real(const int& row, const int& col = 0);
    int         integer(const int& row, const int& col = 0);
    short*      data(void) { return m_data; }
    void        nulval(const short* value);
    short*      nulval(void) { return m_nulval; }

private:
    // Private methods
    void                init_members(void);
    void                copy_members(const GFitsTableShortCol& column);
    void                free_members(void);
    GFitsTableShortCol* clone(void) const;
    std::string         ascii_format(void) const;
    std::string         binary_format(void) const;
    void                alloc_data(void);
    void                release_data(void);
    void                alloc_nulval(const short* value);
    void                init_data(void);
    void*               ptr_data(void) { return m_data; }
    void*               ptr_nulval(void) { return m_nulval; }

    // Private data area
    short* m_data;       //!< Data vector
    short* m_nulval;     //!< NULL value
};

#endif /* GFITSTABLESHORTCOL_HPP */
