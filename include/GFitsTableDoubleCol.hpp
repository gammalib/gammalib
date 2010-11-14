/***************************************************************************
 *        GFitsTableDoubleCol.hpp  - FITS table double column class        *
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
 * @file GFitsTableDoubleCol.hpp
 * @brief GFitsTableDoubleCol class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSTABLEDOUBLECOL_HPP
#define GFITSTABLEDOUBLECOL_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsTableCol.hpp"


/***********************************************************************//**
 * @class GFitsTableDoubleCol
 *
 * @brief Interface for FITS table double column
 *
 * This class implements a FITS table double column.
 ***************************************************************************/
class GFitsTableDoubleCol : public GFitsTableCol {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GFitsTableDoubleCol& column);

public:
    // Constructors and destructors
    GFitsTableDoubleCol(void);
    GFitsTableDoubleCol(const std::string& name, const int& length,
                        const int& size = 1);
    GFitsTableDoubleCol(const GFitsTableDoubleCol& column);
    virtual ~GFitsTableDoubleCol(void);

    // Operators
    GFitsTableDoubleCol& operator= (const GFitsTableDoubleCol& column);
    double&              operator() (const int& row, const int& inx = 0);
    const double&        operator() (const int& row, const int& inx = 0) const;

    // Methods
    std::string string(const int& row, const int& col = 0);
    double      real(const int& row, const int& col = 0);
    int         integer(const int& row, const int& col = 0);
    double*     data(void) { return m_data; }
    void        nullval(const double* value);
    double*     nullval(void) { return m_nulval; }

private:
    // Private methods
    void                 init_members(void);
    void                 copy_members(const GFitsTableDoubleCol& column);
    void                 free_members(void);
    GFitsTableDoubleCol* clone(void) const;
    std::string          ascii_format(void) const;
    std::string          binary_format(void) const;
    void                 alloc_data(void);
    void                 release_data(void);
    void                 alloc_nulval(const double* value);
    void                 init_data(void);
    void*                ptr_data(void) { return m_data; }
    void*                ptr_nulval(void) { return m_nulval; }

    // Private data area
    double* m_data;       //!< Data vector
    double* m_nulval;     //!< NULL value
};

#endif /* GFITSTABLEDOUBLECOL_HPP */
