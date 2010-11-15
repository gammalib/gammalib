/***************************************************************************
 *              GFitsImage.hpp  - FITS image abstract base class           *
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
 * @file GFitsImage.hpp
 * @brief GFitsImage class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSIMAGE_HPP
#define GFITSIMAGE_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsData.hpp"


/***********************************************************************//**
 * @class GFitsImage
 *
 * @brief Abstract interface for the FITS image classes.
 *
 * Implements an abstract interface for FITS images.
 ***************************************************************************/
class GFitsImage : public GFitsData {

    // Friend classes
    friend class GFitsHDU;

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GFitsImage& image);

public:
    // Constructors and destructors
    GFitsImage(void);
    GFitsImage(int naxis, const int* naxes);
    GFitsImage(const GFitsImage& image);
    virtual ~GFitsImage(void);

    // Operators
    GFitsImage& operator= (const GFitsImage& image);

    // Methods
    int bitpix(void) const;
    int naxis(void) const;
    int naxes(int axis) const;
    int num_pixels(void) const;

protected:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsImage& image);
    void free_members(void);
    void open_image(void* vptr);
    void load_image(int datatype, const void* pixels, const void* nulval,
                    int* anynul);
    void save_image(int datatype, const void* pixels);

    // Pure virtual methods
    virtual void        open(void* vptr) = 0;
    virtual void        save(void) = 0;
    virtual void        close(void) = 0;
    virtual GFitsImage* clone(void) const = 0;

    // Private data area
    int   m_bitpix;      //!< Number of Bits/pixel
    int   m_naxis;       //!< Image dimension
    long* m_naxes;       //!< Number of pixels in each dimension
    int   m_num_pixels;  //!< Number of image pixels
};

#endif /* GFITSIMAGE_HPP */
