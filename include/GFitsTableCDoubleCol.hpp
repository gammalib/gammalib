/***************************************************************************
 *  GFitsTableCDoubleCol.hpp  - FITS table double precision complex column *
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
 * @file GFitsTableCDoubleCol.hpp
 * @brief GFitsTableCDoubleCol class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSTABLECDOUBLECOL_HPP
#define GFITSTABLECDOUBLECOL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GFits.hpp"
#include "GFitsTableCol.hpp"


/***********************************************************************//**
 * @class GFitsTableCDoubleCol
 *
 * @brief Interface for FITS table float column
 *
 * This class implements a FITS table float column.
 ***************************************************************************/
class GFitsTableCDoubleCol : public GFitsTableCol {

public:
    // Constructors and destructors
    GFitsTableCDoubleCol(void);
    GFitsTableCDoubleCol(const std::string& name, const int& length,
                         const int& size = 1);
    GFitsTableCDoubleCol(const GFitsTableCDoubleCol& column);
    virtual ~GFitsTableCDoubleCol(void);

    // Operators
    GFitsTableCDoubleCol& operator= (const GFitsTableCDoubleCol& column);
    GFits::cdouble&       operator() (const int& row, const int& inx = 0);
    const GFits::cdouble& operator() (const int& row, const int& inx = 0) const;

    // Methods
    std::string     string(const int& row, const int& col = 0);
    double          real(const int& row, const int& col = 0);
    int             integer(const int& row, const int& col = 0);
    GFits::cdouble* data(void) { return m_data; }
    void            nulval(const GFits::cdouble* value);
    GFits::cdouble* nulval(void) { return m_nulval; }

private:
    // Private methods
    void                  init_members(void);
    void                  copy_members(const GFitsTableCDoubleCol& column);
    void                  free_members(void);
    GFitsTableCDoubleCol* clone(void) const;
    std::string           ascii_format(void) const;
    std::string           binary_format(void) const;
    void                  alloc_data(void);
    void                  release_data(void);
    void                  alloc_nulval(const GFits::cdouble* value);
    void                  init_data(void);
    void*                 ptr_data(void) { return m_data; }
    void*                 ptr_nulval(void) { return m_nulval; }

    // Private data area
    GFits::cdouble* m_data;       //!< Data vector
    GFits::cdouble* m_nulval;     //!< NULL value
};

#endif /* GFITSTABLECDOUBLECOL_HPP */
