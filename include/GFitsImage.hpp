/***************************************************************************
 *                    GFitsImage.hpp  - FITS image class                   *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsImage.hpp
 * @brief GFitsImage class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSIMAGE_HPP
#define GFITSIMAGE_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsCfitsio.hpp"
#include "GFitsData.hpp"

/* __ Namespaces _________________________________________________________ */


/***********************************************************************//**
 * @class GFitsImage
 *
 * @brief Implements a FITS image
 *
 ***************************************************************************/
class GFitsImage : public GFitsData {

    // I/O friends
    friend ostream& operator<< (ostream& os, const GFitsImage& image);

public:
    // Constructors and destructors
    GFitsImage();
    GFitsImage(int naxis, const int* naxes, int bitpix = -64);
    GFitsImage(const GFitsImage& image);
    ~GFitsImage();

    // Operators
    GFitsImage& operator= (const GFitsImage& image);

    // Methods
    void        open(__fitsfile* fptr);
    void        save(void);
    void        close(void);
    void        link(const double* pixels);
    GFitsImage* clone(void) const;
    
private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsImage& image);
    void free_members(void);
    void connect(__fitsfile* fptr);

    // Private data area
    __fitsfile m_fitsfile;    // FITS file
    int        m_bitpix;      // Number of Bits/pixel
    int        m_naxis;       // Image dimension
    long*      m_naxes;       // Number of pixels in each dimension
    int        m_linked;      // Pixels are linked (don't delete them!)
    int        m_type;        // Pixel type
    void*      m_pixels;      // Pixels
};


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/
inline
void GFitsImage::link(const double* pixels)
{
    m_linked = 1;
    m_type   = __TDOUBLE;
    m_pixels = (double*)pixels;
}
inline 
GFitsImage* GFitsImage::clone(void) const 
{
    return new GFitsImage(*this);
}

#endif /* GFITSIMAGE_HPP */
