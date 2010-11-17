/***************************************************************************
 *         GFitsImageFlt.hpp  - FITS single precision image class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsImageFlt.hpp
 * @brief GFitsImageFlt class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSIMAGEFLT_HPP
#define GFITSIMAGEFLT_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsImage.hpp"


/***********************************************************************//**
 * @class GFitsImageFlt
 *
 * @brief Implements a FITS double precision image
 *
 ***************************************************************************/
class GFitsImageFlt : public GFitsImage {

public:
    // Constructors and destructors
    GFitsImageFlt(void);
    GFitsImageFlt(int naxis, const int* naxes);
    GFitsImageFlt(int naxis, const int* naxes, const float* pixels);
    GFitsImageFlt(const GFitsImageFlt& image);
    virtual ~GFitsImageFlt(void);

    // Operators
    GFitsImageFlt& operator= (const GFitsImageFlt& image);
    float&         operator() (const int& ix);
    const float&   operator() (const int& ix) const;
    float&         operator() (const int& ix, const int& iy);
    const float&   operator() (const int& ix, const int& iy) const;
    float&         operator() (const int& ix, const int& iy, const int& iz);
    const float&   operator() (const int& ix, const int& iy, const int& iz) const;
    float&         operator() (const int& ix, const int& iy, const int& iz, const int& it);
    const float&   operator() (const int& ix, const int& iy, const int& iz, const int& it) const;

    // Methods
    void           link(float* pixels);
    void           nulval(const float* value);
    void*          pixels(void);
    double         pixel(const int& ix) { return double((*(this))(ix)); }
    double         pixel(const int& ix, const int& iy) { return double((*(this))(ix, iy)); }
    GFitsImageFlt* clone(void) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsImageFlt& image);
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
    int    m_linked;        // Pixels are linked (don't delete them!)
    float* m_pixels;        // Pixels
    float* m_nulval;        // NULL value
    int    m_anynul;        // Number of NULLs encountered
};

#endif /* GFITSIMAGEFLT_HPP */
