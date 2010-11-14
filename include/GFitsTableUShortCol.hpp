/***************************************************************************
 *    GFitsTableUShortCol.hpp  - FITS table unsigned short column class    *
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
 * @file GFitsTableUShortCol.hpp
 * @brief GFitsTableUShortCol class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSTABLEUSHORTCOL_HPP
#define GFITSTABLEUSHORTCOL_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsTableCol.hpp"


/***********************************************************************//**
 * @class GFitsTableUShortCol
 *
 * @brief Interface for FITS table unsigned short integer column
 *
 * This class implements a FITS table unsigned short integer column.
 ***************************************************************************/
class GFitsTableUShortCol : public GFitsTableCol {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GFitsTableUShortCol& column);

public:
    // Constructors and destructors
    GFitsTableUShortCol(void);
    GFitsTableUShortCol(const std::string& name, const int& length,
                        const int& size = 1);
    GFitsTableUShortCol(const GFitsTableUShortCol& column);
    virtual ~GFitsTableUShortCol(void);

    // Operators
    GFitsTableUShortCol&  operator= (const GFitsTableUShortCol& column);
    unsigned short&       operator() (const int& row, const int& inx = 0);
    const unsigned short& operator() (const int& row, const int& inx = 0) const;

    // Methods
    std::string     string(const int& row, const int& col = 0);
    double          real(const int& row, const int& col = 0);
    int             integer(const int& row, const int& col = 0);
    unsigned short* data(void) { return m_data; }
    void            nullval(const unsigned short* value);
    unsigned short* nullval(void) { return m_nulval; }

private:
    // Private methods
    void                 init_members(void);
    void                 copy_members(const GFitsTableUShortCol& column);
    void                 free_members(void);
    GFitsTableUShortCol* clone(void) const;
    std::string          ascii_format(void) const;
    std::string          binary_format(void) const;
    void                 alloc_data(void);
    void                 release_data(void);
    void                 alloc_nulval(const unsigned short* value);
    void                 init_data(void);
    void*                ptr_data(void) { return m_data; }
    void*                ptr_nulval(void) { return m_nulval; }

    // Private data area
    unsigned short* m_data;       //!< Data vector
    unsigned short* m_nulval;     //!< NULL value
};

#endif /* GFITSTABLEUSHORTCOL_HPP */
