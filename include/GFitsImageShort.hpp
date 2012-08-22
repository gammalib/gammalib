/***************************************************************************
 *          GFitsImageShort.hpp  - FITS short integer image class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
 * @file GFitsImageShort.hpp
 * @brief GFitsImageShort class definition.
 * @author J. Knoedlseder
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
    explicit GFitsImageShort(int nx, const short* pixels = NULL);
    explicit GFitsImageShort(int nx, int ny, const short* pixels = NULL);
    explicit GFitsImageShort(int nx, int ny, int nz, const short* pixels = NULL);
    explicit GFitsImageShort(int nx, int ny, int nz, int nt, const short* pixels = NULL);
    explicit GFitsImageShort(int naxis, const int* naxes, const short* pixels = NULL);
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
    int              type(void) const;
    GFitsImageShort* clone(void) const;

private:
    // Private methods
    void  init_members(void);
    void  copy_members(const GFitsImageShort& image);
    void  free_members(void);
    void  alloc_data(void);
    void  init_data(void);
    void  release_data(void);
    void  construct_data(const short* pixels);
    void  load_data(void) const;
    void  alloc_nulval(const void* value);
    void* ptr_data(void) { return m_pixels; }
    void* ptr_nulval(void) { return m_nulval; }

    // Private data area
    short* m_pixels;      //!< Pixels
    short* m_nulval;      //!< NULL value
};

#endif /* GFITSIMAGESHORT_HPP */
