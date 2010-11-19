/***************************************************************************
 *      GFitsImageLongLong.hpp  - FITS long long integer image class       *
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
 * @file GFitsImageLongLong.hpp
 * @brief GFitsImageLongLong class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSIMAGELONGLONG_HPP
#define GFITSIMAGELONGLONG_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsImage.hpp"


/***********************************************************************//**
 * @class GFitsImageLongLong
 *
 * @brief Implements a FITS long long integer image
 ***************************************************************************/
class GFitsImageLongLong : public GFitsImage {

public:
    // Constructors and destructors
    GFitsImageLongLong(void);
    explicit GFitsImageLongLong(int nx, const long long* pixels = NULL);
    explicit GFitsImageLongLong(int nx, int ny, const long long* pixels = NULL);
    explicit GFitsImageLongLong(int nx, int ny, int nz, const long long* pixels = NULL);
    explicit GFitsImageLongLong(int nx, int ny, int nz, int nt, const long long* pixels = NULL);
    explicit GFitsImageLongLong(int naxis, const int* naxes, const long long* pixels = NULL);
    GFitsImageLongLong(const GFitsImageLongLong& image);
    virtual ~GFitsImageLongLong(void);

    // Operators
    GFitsImageLongLong& operator= (const GFitsImageLongLong& image);
    long long&          operator() (const int& ix);
    long long&          operator() (const int& ix, const int& iy);
    long long&          operator() (const int& ix, const int& iy, const int& iz);
    long long&          operator() (const int& ix, const int& iy, const int& iz, const int& it);
    const long long&    operator() (const int& ix) const;
    const long long&    operator() (const int& ix, const int& iy) const;
    const long long&    operator() (const int& ix, const int& iy, const int& iz) const;
    const long long&    operator() (const int& ix, const int& iy, const int& iz, const int& it) const;

    // Methods
    long long&          at(const int& ix);
    long long&          at(const int& ix, const int& iy);
    long long&          at(const int& ix, const int& iy, const int& iz);
    long long&          at(const int& ix, const int& iy, const int& iz, const int& it);
    const long long&    at(const int& ix) const;
    const long long&    at(const int& ix, const int& iy) const;
    const long long&    at(const int& ix, const int& iy, const int& iz) const;
    const long long&    at(const int& ix, const int& iy, const int& iz, const int& it) const;
    double              pixel(const int& ix) const;
    double              pixel(const int& ix, const int& iy) const;
    double              pixel(const int& ix, const int& iy, const int& iz) const;
    double              pixel(const int& ix, const int& iy, const int& iz, const int& it) const;
    void*               pixels(void);
    GFitsImageLongLong* clone(void) const;

private:
    // Private methods
    void  init_members(void);
    void  copy_members(const GFitsImageLongLong& image);
    void  free_members(void);
    void  alloc_data(void);
    void  init_data(void);
    void  release_data(void);
    void  construct_data(const long long* pixels);
    void  alloc_nulval(const void* value);
    void* ptr_data(void) { return m_pixels; }
    void* ptr_nulval(void) { return m_nulval; }
    int   type(void) const;

    // Private data area
    long long* m_pixels;      //!< Pixels
    long long* m_nulval;      //!< NULL value
};

#endif /* GFITSIMAGELONGLONG_HPP */
