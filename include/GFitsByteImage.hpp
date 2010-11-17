/***************************************************************************
 *               GFitsByteImage.hpp  - FITS Byte image class               *
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
 * @file GFitsByteImage.hpp
 * @brief GFitsByteImage class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSBYTEIMAGE_HPP
#define GFITSBYTEIMAGE_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsImage.hpp"


/***********************************************************************//**
 * @class GFitsByteImage
 *
 * @brief Implements a FITS byte image
 ***************************************************************************/
class GFitsByteImage : public GFitsImage {

public:
    // Constructors and destructors
    GFitsByteImage(void);
    GFitsByteImage(int naxis, const int* naxes);
    GFitsByteImage(const GFitsByteImage& image);
    virtual ~GFitsByteImage(void);

    // Operators
    GFitsByteImage&      operator= (const GFitsByteImage& image);
    unsigned char&       operator() (const int& ix);
    const unsigned char& operator() (const int& ix) const;
    unsigned char&       operator() (const int& ix, const int& iy);
    const unsigned char& operator() (const int& ix, const int& iy) const;
    unsigned char&       operator() (const int& ix, const int& iy, const int& iz);
    const unsigned char& operator() (const int& ix, const int& iy, const int& iz) const;
    unsigned char&       operator() (const int& ix, const int& iy, const int& iz, const int& it);
    const unsigned char& operator() (const int& ix, const int& iy, const int& iz, const int& it) const;

    // Methods
    void            nulval(const unsigned char* value);
    void*           pixels(void);
    double          pixel(const int& ix) const { return double((*(this))(ix)); }
    double          pixel(const int& ix, const int& iy) const { return double((*(this))(ix, iy)); }
    double          pixel(const int& ix, const int& iy, const int& iz) const { return double((*(this))(ix, iy, iz)); }
    double          pixel(const int& ix, const int& iy, const int& iz, const int& it) const { return double((*(this))(ix, iy, iz, it)); }
    GFitsByteImage* clone(void) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsByteImage& image);
    void free_members(void);
    void fetch_pixels(void);
    int  type(void) const; 

    // Stuff to allow for compilation (new GFitsImage interface)
    void  alloc_data(void);
    void  init_data(void);
    void  release_data(void);
    void  alloc_nulval(const void* value);
    void* ptr_data(void) { return m_pixels; }
    void* ptr_nulval(void) { return m_nulval; }

    // Private data area
    unsigned char* m_pixels;      //!< Image pixels
    unsigned char* m_nulval;      //!< NULL value
    int            m_anynul;      //!< Number of NULLs encountered
};

#endif /* GFITSBYTEIMAGE_HPP */
