/***************************************************************************
 *         GFitsImageFloat.hpp - Single precision FITS image class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsImageFloat.hpp
 * @brief Single precision FITS image class definition
 * @author Juergen Knoedlseder
 */

#ifndef GFITSIMAGEFLOAT_HPP
#define GFITSIMAGEFLOAT_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsImage.hpp"


/***********************************************************************//**
 * @class GFitsImageFloat
 *
 * @brief Single precision FITS image class
 ***************************************************************************/
class GFitsImageFloat : public GFitsImage {

public:
    // Constructors and destructors
    GFitsImageFloat(void);
    GFitsImageFloat(const int& nx, const float* pixels = NULL);
    GFitsImageFloat(const int& nx, const int& ny, const float* pixels = NULL);
    GFitsImageFloat(const int& nx, const int& ny, const int& nz, const float* pixels = NULL);
    GFitsImageFloat(const int& nx, const int& ny, const int& nz, const int& nt, const float* pixels = NULL);
    GFitsImageFloat(const int& naxis, const int* naxes, const float* pixels = NULL);
    GFitsImageFloat(const GFitsImageFloat& image);
    virtual ~GFitsImageFloat(void);

    // Operators
    GFitsImageFloat& operator=(const GFitsImageFloat& image);
    float&           operator()(const int& ix);
    float&           operator()(const int& ix, const int& iy);
    float&           operator()(const int& ix, const int& iy, const int& iz);
    float&           operator()(const int& ix, const int& iy, const int& iz, const int& it);
    const float&     operator()(const int& ix) const;
    const float&     operator()(const int& ix, const int& iy) const;
    const float&     operator()(const int& ix, const int& iy, const int& iz) const;
    const float&     operator()(const int& ix, const int& iy, const int& iz, const int& it) const;

    // Methods
    void             clear(void);
    GFitsImageFloat* clone(void) const;
    float&           at(const int& ix);
    float&           at(const int& ix, const int& iy);
    float&           at(const int& ix, const int& iy, const int& iz);
    float&           at(const int& ix, const int& iy, const int& iz, const int& it);
    const float&     at(const int& ix) const;
    const float&     at(const int& ix, const int& iy) const;
    const float&     at(const int& ix, const int& iy, const int& iz) const;
    const float&     at(const int& ix, const int& iy, const int& iz, const int& it) const;
    double           pixel(const int& ix) const;
    double           pixel(const int& ix, const int& iy) const;
    double           pixel(const int& ix, const int& iy, const int& iz) const;
    double           pixel(const int& ix, const int& iy, const int& iz, const int& it) const;
    void*            pixels(void);
    int              type(void) const;

private:
    // Private methods
    void  init_members(void);
    void  copy_members(const GFitsImageFloat& image);
    void  free_members(void);
    void  alloc_data(void);
    void  init_data(void);
    void  release_data(void);
    void  construct_data(const float* pixels);
    void  load_data(void) const;
    void  alloc_nulval(const void* value);
    void* ptr_data(void) { return m_pixels; }
    void* ptr_nulval(void) { return m_nulval; }

    // Private data area
    float* m_pixels;      //!< Pixels
    float* m_nulval;      //!< NULL value
};

#endif /* GFITSIMAGEFLOAT_HPP */
