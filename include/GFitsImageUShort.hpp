/***************************************************************************
 *         GFitsImageUShort.hpp  - FITS unsigned short image class         *
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
 * @file GFitsImageUShort.hpp
 * @brief GFitsImageUShort class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSIMAGEUSHORT_HPP
#define GFITSIMAGEUSHORT_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsImage.hpp"


/***********************************************************************//**
 * @class GFitsImageUShort
 *
 * @brief Implements a FITS unsigned short integer image
 ***************************************************************************/
class GFitsImageUShort : public GFitsImage {

public:
    // Constructors and destructors
    GFitsImageUShort(void);
    GFitsImageUShort(int naxis, const int* naxes, const unsigned short* pixels = NULL);
    GFitsImageUShort(const GFitsImageUShort& image);
    virtual ~GFitsImageUShort(void);

    // Operators
    GFitsImageUShort&     operator= (const GFitsImageUShort& image);
    unsigned short&       operator() (const int& ix);
    unsigned short&       operator() (const int& ix, const int& iy);
    unsigned short&       operator() (const int& ix, const int& iy, const int& iz);
    unsigned short&       operator() (const int& ix, const int& iy, const int& iz, const int& it);
    const unsigned short& operator() (const int& ix) const;
    const unsigned short& operator() (const int& ix, const int& iy) const;
    const unsigned short& operator() (const int& ix, const int& iy, const int& iz) const;
    const unsigned short& operator() (const int& ix, const int& iy, const int& iz, const int& it) const;

    // Methods
    unsigned short&       at(const int& ix);
    unsigned short&       at(const int& ix, const int& iy);
    unsigned short&       at(const int& ix, const int& iy, const int& iz);
    unsigned short&       at(const int& ix, const int& iy, const int& iz, const int& it);
    const unsigned short& at(const int& ix) const;
    const unsigned short& at(const int& ix, const int& iy) const;
    const unsigned short& at(const int& ix, const int& iy, const int& iz) const;
    const unsigned short& at(const int& ix, const int& iy, const int& iz, const int& it) const;
    double                pixel(const int& ix) const;
    double                pixel(const int& ix, const int& iy) const;
    double                pixel(const int& ix, const int& iy, const int& iz) const;
    double                pixel(const int& ix, const int& iy, const int& iz, const int& it) const;
    void*                 pixels(void);
    GFitsImageUShort*     clone(void) const;

private:
    // Private methods
    void  init_members(void);
    void  copy_members(const GFitsImageUShort& image);
    void  free_members(void);
    void  alloc_data(void);
    void  init_data(void);
    void  release_data(void);
    void  alloc_nulval(const void* value);
    void* ptr_data(void) { return m_pixels; }
    void* ptr_nulval(void) { return m_nulval; }
    int   type(void) const;

    // Private data area
    unsigned short* m_pixels;      //!< Pixels
    unsigned short* m_nulval;      //!< NULL value
};

#endif /* GFITSIMAGEUSHORT_HPP */
