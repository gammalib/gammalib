/***************************************************************************
 *        GFitsImageDouble.hpp  - FITS double precision image class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsImageDouble.hpp
 * @brief GFitsImageDouble class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSIMAGEDOUBLE_HPP
#define GFITSIMAGEDOUBLE_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsImage.hpp"


/***********************************************************************//**
 * @class GFitsImageDouble
 *
 * @brief Implements a FITS double precision image
 ***************************************************************************/
class GFitsImageDouble : public GFitsImage {

public:
    // Constructors and destructors
    GFitsImageDouble(void);
    GFitsImageDouble(int naxis, const int* naxes, const double* pixels = NULL);
    GFitsImageDouble(const GFitsImageDouble& image);
    virtual ~GFitsImageDouble(void);

    // Operators
    GFitsImageDouble& operator= (const GFitsImageDouble& image);
    double&           operator() (const int& ix);
    const double&     operator() (const int& ix) const;
    double&           operator() (const int& ix, const int& iy);
    const double&     operator() (const int& ix, const int& iy) const;
    double&           operator() (const int& ix, const int& iy, const int& iz);
    const double&     operator() (const int& ix, const int& iy, const int& iz) const;
    double&           operator() (const int& ix, const int& iy, const int& iz, const int& it);
    const double&     operator() (const int& ix, const int& iy, const int& iz, const int& it) const;

    // Methods
    double&           at(const int& ix);
    const double&     at(const int& ix) const;

    void*             pixels(void);
    double            pixel(const int& ix) { return double((*(this))(ix)); }
    double            pixel(const int& ix, const int& iy) { return double((*(this))(ix, iy)); }
    GFitsImageDouble* clone(void) const;

private:
    // Private methods
    void  init_members(void);
    void  copy_members(const GFitsImageDouble& image);
    void  free_members(void);
    void  alloc_data(void);
    void  init_data(void);
    void  release_data(void);
    void  alloc_nulval(const void* value);
    void* ptr_data(void) { return m_pixels; }
    void* ptr_nulval(void) { return m_nulval; }
    int   type(void) const { return __TDOUBLE; }

    // Private data area
    double* m_pixels;      //!< Pixels
    double* m_nulval;      //!< NULL value
};

#endif /* GFITSIMAGEDOUBLE_HPP */
