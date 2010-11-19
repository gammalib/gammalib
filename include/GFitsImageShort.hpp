/***************************************************************************
 *          GFitsImageShort.hpp  - FITS short integer image class          *
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
 * @file GFitsImageShort.hpp
 * @brief GFitsImageShort class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSIMAGESHORT_HPP
#define GFITSIMAGESHORT_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsImage.hpp"


/***********************************************************************//**
 * @class GFitsImageShort
 *
 * @brief Implements a FITS short integer image
 ***************************************************************************/
class GFitsImageShort : public GFitsImage {

public:
    // Constructors and destructors
    GFitsImageShort(void);
    GFitsImageShort(int naxis, const int* naxes, const short* pixels = NULL);
    GFitsImageShort(const GFitsImageShort& image);
    virtual ~GFitsImageShort(void);

    // Operators
    GFitsImageShort& operator= (const GFitsImageShort& image);
    short&           operator() (const int& ix);
    short&           operator() (const int& ix, const int& iy);
    short&           operator() (const int& ix, const int& iy, const int& iz);
    short&           operator() (const int& ix, const int& iy, const int& iz, const int& it);
    const short&     operator() (const int& ix) const;
    const short&     operator() (const int& ix, const int& iy) const;
    const short&     operator() (const int& ix, const int& iy, const int& iz) const;
    const short&     operator() (const int& ix, const int& iy, const int& iz, const int& it) const;

    // Methods
    short&           at(const int& ix);
    short&           at(const int& ix, const int& iy);
    short&           at(const int& ix, const int& iy, const int& iz);
    short&           at(const int& ix, const int& iy, const int& iz, const int& it);
    const short&     at(const int& ix) const;
    const short&     at(const int& ix, const int& iy) const;
    const short&     at(const int& ix, const int& iy, const int& iz) const;
    const short&     at(const int& ix, const int& iy, const int& iz, const int& it) const;
    double           pixel(const int& ix) const;
    double           pixel(const int& ix, const int& iy) const;
    double           pixel(const int& ix, const int& iy, const int& iz) const;
    double           pixel(const int& ix, const int& iy, const int& iz, const int& it) const;
    void*            pixels(void);
    GFitsImageShort* clone(void) const;

private:
    // Private methods
    void  init_members(void);
    void  copy_members(const GFitsImageShort& image);
    void  free_members(void);
    void  alloc_data(void);
    void  init_data(void);
    void  release_data(void);
    void  alloc_nulval(const void* value);
    void* ptr_data(void) { return m_pixels; }
    void* ptr_nulval(void) { return m_nulval; }
    int   type(void) const;

    // Private data area
    short* m_pixels;      //!< Pixels
    short* m_nulval;      //!< NULL value
};

#endif /* GFITSIMAGESHORT_HPP */
