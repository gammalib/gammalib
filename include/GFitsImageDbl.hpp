/***************************************************************************
 *         GFitsImageDbl.hpp  - FITS double precision image class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsImageDbl.hpp
 * @brief GFitsImageDbl class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSIMAGEDBL_HPP
#define GFITSIMAGEDBL_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsImage.hpp"


/***********************************************************************//**
 * @class GFitsImageDbl
 *
 * @brief Implements a FITS double precision image
 ***************************************************************************/
class GFitsImageDbl : public GFitsImage {

public:
    // Constructors and destructors
    GFitsImageDbl(void);
    GFitsImageDbl(int naxis, const int* naxes);
    GFitsImageDbl(int naxis, const int* naxes, const double* pixels);
    GFitsImageDbl(const GFitsImageDbl& image);
    virtual ~GFitsImageDbl(void);

    // Operators
    GFitsImageDbl& operator= (const GFitsImageDbl& image);
    double&        operator() (const int& ix);
    const double&  operator() (const int& ix) const;
    double&        operator() (const int& ix, const int& iy);
    const double&  operator() (const int& ix, const int& iy) const;
    double&        operator() (const int& ix, const int& iy, const int& iz);
    const double&  operator() (const int& ix, const int& iy, const int& iz)
                               const;
    double&        operator() (const int& ix, const int& iy, const int& iz,
                               const int& it);
    const double&  operator() (const int& ix, const int& iy, const int& iz,
                               const int& it) const;

    // Methods
    void           link(double* pixels);
    void           nulval(const double* value);
    void*          pixels(void);
    double         pixel(const int& ix) { return double((*(this))(ix)); }
    double         pixel(const int& ix, const int& iy) { return double((*(this))(ix, iy)); }
    GFitsImageDbl* clone(void) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsImageDbl& image);
    void free_members(void);
    void fetch_pixels(void);
    int  type(void) const;

    // Stuff to allow for compilation (new GFitsImage interface)
    void  alloc_data(void) { return; }
    void  init_data(void) { return; }
    void  release_data(void) { return; }
    void  alloc_nulval(const void* value) { return; }
    void* ptr_data(void) { void* p; return p; }
    void* ptr_nulval(void) { void* p; return p; }

    // Private data area
    int     m_linked;        // Pixels are linked (don't delete them!)
    double* m_pixels;        // Pixels
    double* m_nulval;        // NULL value
    int     m_anynul;        // Number of NULLs encountered
};

#endif /* GFITSIMAGEDBL_HPP */
