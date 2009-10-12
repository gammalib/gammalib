/***************************************************************************
 *         GFitsImageFlt.hpp  - FITS single precision image class          *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2009 by Jurgen Knodlseder                   *
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
#include "GFitsCfitsio.hpp"
#include "GFitsData.hpp"
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
    GFitsImageFlt();
    GFitsImageFlt(int naxis, const int* naxes);
    GFitsImageFlt(int naxis, const int* naxes, const float* pixels);
    GFitsImageFlt(const GFitsImageFlt& image);
    ~GFitsImageFlt();

    // Operators
    GFitsImageFlt& operator= (const GFitsImageFlt& image);
    float&         operator() (const int& ix);
    const float&   operator() (const int& ix) const;
    float&         operator() (const int& ix, const int& iy);
    const float&   operator() (const int& ix, const int& iy) const;
    float&         operator() (const int& ix, const int& iy, const int& iz);
    const float&   operator() (const int& ix, const int& iy, const int& iz)
                               const;
    float&         operator() (const int& ix, const int& iy, const int& iz,
                               const int& it);
    const float&   operator() (const int& ix, const int& iy, const int& iz,
                               const int& it) const;

    // Methods
    void   link(float* pixels);
    void   set_nullval(const float* value);
    float* pixels(void);

private:
    // Private methods
    void           init_members(void);
    void           copy_members(const GFitsImageFlt& image);
    void           free_members(void);
    void           fetch_pixels(void);
    void           open(__fitsfile* fptr);
    void           save(void);
    void           close(void);
    GFitsImageFlt* clone(void) const;

    // Private data area
    int    m_linked;        // Pixels are linked (don't delete them!)
    float* m_pixels;        // Pixels
    float* m_nulval;        // NULL value
    int    m_anynul;        // Number of NULLs encountered
};

#endif /* GFITSIMAGEFLT_HPP */
