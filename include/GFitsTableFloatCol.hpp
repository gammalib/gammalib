/***************************************************************************
 *         GFitsTableFloatCol.hpp  - FITS table float column class         *
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
 * @file GFitsTableFloatCol.hpp
 * @brief GFitsTableFloatCol class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSTABLEFLOATCOL_HPP
#define GFITSTABLEFLOATCOL_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsTableCol.hpp"


/***********************************************************************//**
 * @class GFitsTableFloatCol
 *
 * @brief Interface for FITS table float column
 *
 * This class implements a FITS table float column.
 ***************************************************************************/
class GFitsTableFloatCol : public GFitsTableCol {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GFitsTableFloatCol& column);

public:
    // Constructors and destructors
    GFitsTableFloatCol(void);
    GFitsTableFloatCol(const std::string& name, const int& length,
                       const int& size = 1);
    GFitsTableFloatCol(const GFitsTableFloatCol& column);
    virtual ~GFitsTableFloatCol(void);

    // Operators
    GFitsTableFloatCol& operator= (const GFitsTableFloatCol& column);
    float&              operator() (const int& row, const int& inx = 0);
    const float&        operator() (const int& row, const int& inx = 0) const;

    // Methods
    std::string string(const int& row, const int& col = 0);
    double      real(const int& row, const int& col = 0);
    int         integer(const int& row, const int& col = 0);
    float*      data(void) { return m_data; }
    void        nullval(const float* value);
    float*      nullval(void) { return m_nulval; }

private:
    // Private methods
    void                init_members(void);
    void                copy_members(const GFitsTableFloatCol& column);
    void                free_members(void);
    GFitsTableFloatCol* clone(void) const;
    std::string         ascii_format(void) const;
    std::string         binary_format(void) const;
    void                alloc_data(void);
    void                release_data(void);
    void                alloc_nulval(const float* value);
    void                init_data(void);
    void*               ptr_data(void) { return m_data; }
    void*               ptr_nulval(void) { return m_nulval; }

    // Private data area
    float* m_data;       //!< Data vector
    float* m_nulval;     //!< NULL value
};

#endif /* GFITSTABLEFLOATCOL_HPP */
