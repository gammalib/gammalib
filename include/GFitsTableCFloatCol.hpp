/***************************************************************************
 *  GFitsTableCFloatCol.hpp  - FITS table single precision complex column  *
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
 * @file GFitsTableCFloatCol.hpp
 * @brief GFitsTableCFloatCol class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSTABLECFLOATCOL_HPP
#define GFITSTABLECFLOATCOL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GFits.hpp"
#include "GFitsTableCol.hpp"


/***********************************************************************//**
 * @class GFitsTableCFloatCol
 *
 * @brief Interface for FITS table float column
 *
 * This class implements a FITS table float column.
 ***************************************************************************/
class GFitsTableCFloatCol : public GFitsTableCol {

public:
    // Constructors and destructors
    GFitsTableCFloatCol(void);
    GFitsTableCFloatCol(const std::string& name, const int& length,
                         const int& size = 1);
    GFitsTableCFloatCol(const GFitsTableCFloatCol& column);
    virtual ~GFitsTableCFloatCol(void);

    // Operators
    GFitsTableCFloatCol& operator= (const GFitsTableCFloatCol& column);
    GFits::cfloat&        operator() (const int& row, const int& inx = 0);
    const GFits::cfloat&  operator() (const int& row, const int& inx = 0) const;

    // Methods
    std::string    string(const int& row, const int& col = 0);
    double         real(const int& row, const int& col = 0);
    int            integer(const int& row, const int& col = 0);
    GFits::cfloat* data(void) { return m_data; }
    void           nulval(const GFits::cfloat* value);
    GFits::cfloat* nulval(void) { return m_nulval; }

private:
    // Private methods
    void                  init_members(void);
    void                  copy_members(const GFitsTableCFloatCol& column);
    void                  free_members(void);
    GFitsTableCFloatCol* clone(void) const;
    std::string           ascii_format(void) const;
    std::string           binary_format(void) const;
    void                  alloc_data(void);
    void                  release_data(void);
    void                  alloc_nulval(const GFits::cfloat* value);
    void                  init_data(void);
    void*                 ptr_data(void) { return m_data; }
    void*                 ptr_nulval(void) { return m_nulval; }

    // Private data area
    GFits::cfloat* m_data;       //!< Data vector
    GFits::cfloat* m_nulval;     //!< NULL value
};

#endif /* GFITSTABLECFLOATCOL_HPP */
