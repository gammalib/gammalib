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
#include "GFitsHDU.hpp"


/***********************************************************************//**
 * @class GFitsImage
 *
 * @brief Abstract interface for the FITS image classes.
 *
 * This class defines the abstract interface for a FITS image.
 ***************************************************************************/
class GFitsImage : public GFitsHDU {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GFitsImage& image);

public:
    // Constructors and destructors
    GFitsImage(void);
    GFitsImage(int bitpix, int naxis, const int* naxes);
    GFitsImage(const GFitsImage& image);
    virtual ~GFitsImage(void);

    // Operators
    GFitsImage& operator= (const GFitsImage& image);

    // Pure virtual methods
    virtual void*       pixels(void) = 0;
    virtual GFitsImage* clone(void) const = 0;

    // Implemented pure virtual methods
    HDUType exttype(void) const { return HT_IMAGE; }

    // Methods
    virtual int bitpix(void) const;
    virtual int naxis(void) const;
    virtual int naxes(int axis) const;
    virtual int num_pixels(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GFitsImage& image);
    void free_members(void);
    void init_image_header(void);
    void data_open(void* vptr);
    void data_save(void);
    void data_close(void);
    void data_connect(void* vptr);
    void open_image(void* vptr);
    void load_image(int datatype, const void* pixels, const void* nulval,
                    int* anynul);
    void save_image(int datatype, const void* pixels);

    // Pure virtual protected methods
    virtual int type(void) const = 0;

    // Protected data area
    int   m_bitpix;      //!< Number of Bits/pixel
    int   m_naxis;       //!< Image dimension
    long* m_naxes;       //!< Number of pixels in each dimension
    int   m_num_pixels;  //!< Number of image pixels
};

#endif /* GFITSIMAGE_HPP */
