/***************************************************************************
 *        GFitsImageDouble.hpp  - FITS double precision image class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2012 by Juergen Knoedlseder                         *
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
 * @file GFitsImageDouble.hpp
 * @brief GFitsImageDouble class definition.
 * @author J. Knoedlseder
 */

#ifndef GFITSIMAGEDOUBLE_HPP
#define GFITSIMAGEDOUBLE_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsImage.hpp"


/***********************************************************************//**
 * @class GFitsImageDouble
 *
 * @brief Implements a FITS double precision image
 ***************************************************************************/
class GFitsImageDouble : public GFitsImage {

public:
    // Constructors and destructors
    GFitsImageDouble(void);
    explicit GFitsImageDouble(int nx, const double* pixels = NULL);
    explicit GFitsImageDouble(int nx, int ny, const double* pixels = NULL);
    explicit GFitsImageDouble(int nx, int ny, int nz, const double* pixels = NULL);
    explicit GFitsImageDouble(int nx, int ny, int nz, int nt, const double* pixels = NULL);
    explicit GFitsImageDouble(int naxis, const int* naxes, const double* pixels = NULL);
    GFitsImageDouble(const GFitsImageDouble& image);
    virtual ~GFitsImageDouble(void);

    // Operators
    GFitsImageDouble& operator= (const GFitsImageDouble& image);
    double&           operator() (const int& ix);
    double&           operator() (const int& ix, const int& iy);
    double&           operator() (const int& ix, const int& iy, const int& iz);
    double&           operator() (const int& ix, const int& iy, const int& iz, const int& it);
    const double&     operator() (const int& ix) const;
    const double&     operator() (const int& ix, const int& iy) const;
    const double&     operator() (const int& ix, const int& iy, const int& iz) const;
    const double&     operator() (const int& ix, const int& iy, const int& iz, const int& it) const;

    // Methods
    double&           at(const int& ix);
    double&           at(const int& ix, const int& iy);
    double&           at(const int& ix, const int& iy, const int& iz);
    double&           at(const int& ix, const int& iy, const int& iz, const int& it);
    const double&     at(const int& ix) const;
    const double&     at(const int& ix, const int& iy) const;
    const double&     at(const int& ix, const int& iy, const int& iz) const;
    const double&     at(const int& ix, const int& iy, const int& iz, const int& it) const;
    double            pixel(const int& ix) const;
    double            pixel(const int& ix, const int& iy) const;
    double            pixel(const int& ix, const int& iy, const int& iz) const;
    double            pixel(const int& ix, const int& iy, const int& iz, const int& it) const;
    void*             pixels(void);
    int               type(void) const;
    GFitsImageDouble* clone(void) const;

private:
    // Private methods
    void  init_members(void);
    void  copy_members(const GFitsImageDouble& image);
    void  free_members(void);
    void  alloc_data(void);
    void  init_data(void);
    void  release_data(void);
    void  construct_data(const double* pixels);
    void  load_data(void) const;
    void  alloc_nulval(const void* value);
    void* ptr_data(void) { return m_pixels; }
    void* ptr_nulval(void) { return m_nulval; }

    // Private data area
    double* m_pixels;      //!< Pixels
    double* m_nulval;      //!< NULL value
};

#endif /* GFITSIMAGEDOUBLE_HPP */
