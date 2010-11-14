/***************************************************************************
 *     GFitsTableLongLongCol.hpp  - FITS table long long column class      *
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
 * @file GFitsTableLongLongCol.hpp
 * @brief GFitsTableLongLongCol class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSTABLELONGLONGCOL_HPP
#define GFITSTABLELONGLONGCOL_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsTableCol.hpp"


/***********************************************************************//**
 * @class GFitsTableLongLongCol
 *
 * @brief Interface for FITS table long integer column
 *
 * This class implements a FITS table long integer column.
 ***************************************************************************/
class GFitsTableLongLongCol : public GFitsTableCol {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GFitsTableLongLongCol& column);

public:
    // Constructors and destructors
    GFitsTableLongLongCol(void);
    GFitsTableLongLongCol(const std::string& name, const int& length,
                          const int& size = 1);
    GFitsTableLongLongCol(const GFitsTableLongLongCol& column);
    virtual ~GFitsTableLongLongCol(void);

    // Operators
    GFitsTableLongLongCol& operator= (const GFitsTableLongLongCol& column);
    long long&             operator() (const int& row, const int& inx = 0);
    const long long&       operator() (const int& row, const int& inx = 0) const;

    // Methods
    std::string string(const int& row, const int& col = 0);
    double      real(const int& row, const int& col = 0);
    int         integer(const int& row, const int& col = 0);
    long long*  data(void) { return m_data; }
    void        nullval(const long long* value);
    long long*  nullval(void) { return m_nulval; }

private:
    // Private methods
    void                   init_members(void);
    void                   copy_members(const GFitsTableLongLongCol& column);
    void                   free_members(void);
    GFitsTableLongLongCol* clone(void) const;
    std::string            ascii_format(void) const;
    std::string            binary_format(void) const;
    void                   alloc_data(void);
    void                   release_data(void);
    void                   alloc_nulval(const long long* value);
    void                   init_data(void);
    void*                  ptr_data(void) { return m_data; }
    void*                  ptr_nulval(void) { return m_nulval; }

    // Private data area
    long long* m_data;       //!< Data vector
    long long* m_nulval;     //!< NULL value
};

#endif /* GFITSTABLELONGLONGCOL_HPP */
