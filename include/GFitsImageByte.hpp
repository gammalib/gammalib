/***************************************************************************
 *               GFitsImageByte.hpp  - FITS Byte image class               *
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
 * @file GFitsImageByte.hpp
 * @brief GFitsImageByte class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSIMAGEBYTE_HPP
#define GFITSIMAGEBYTE_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsImage.hpp"


/***********************************************************************//**
 * @class GFitsImageByte
 *
 * @brief Implements a FITS Byte image
 ***************************************************************************/
class GFitsImageByte : public GFitsImage {

public:
    // Constructors and destructors
    GFitsImageByte(void);
    explicit GFitsImageByte(int nx, const unsigned char* pixels = NULL);
    explicit GFitsImageByte(int nx, int ny, const unsigned char* pixels = NULL);
    explicit GFitsImageByte(int nx, int ny, int nz, const unsigned char* pixels = NULL);
    explicit GFitsImageByte(int nx, int ny, int nz, int nt, const unsigned char* pixels = NULL);
    explicit GFitsImageByte(int naxis, const int* naxes, const unsigned char* pixels = NULL);
    GFitsImageByte(const GFitsImageByte& image);
    virtual ~GFitsImageByte(void);

    // Operators
    GFitsImageByte&      operator= (const GFitsImageByte& image);
    unsigned char&       operator() (const int& ix);
    unsigned char&       operator() (const int& ix, const int& iy);
    unsigned char&       operator() (const int& ix, const int& iy, const int& iz);
    unsigned char&       operator() (const int& ix, const int& iy, const int& iz, const int& it);
    const unsigned char& operator() (const int& ix) const;
    const unsigned char& operator() (const int& ix, const int& iy) const;
    const unsigned char& operator() (const int& ix, const int& iy, const int& iz) const;
    const unsigned char& operator() (const int& ix, const int& iy, const int& iz, const int& it) const;

    // Methods
    void                 clear(void);
    GFitsImageByte*      clone(void) const;
    unsigned char&       at(const int& ix);
    unsigned char&       at(const int& ix, const int& iy);
    unsigned char&       at(const int& ix, const int& iy, const int& iz);
    unsigned char&       at(const int& ix, const int& iy, const int& iz, const int& it);
    const unsigned char& at(const int& ix) const;
    const unsigned char& at(const int& ix, const int& iy) const;
    const unsigned char& at(const int& ix, const int& iy, const int& iz) const;
    const unsigned char& at(const int& ix, const int& iy, const int& iz, const int& it) const;
    double               pixel(const int& ix) const;
    double               pixel(const int& ix, const int& iy) const;
    double               pixel(const int& ix, const int& iy, const int& iz) const;
    double               pixel(const int& ix, const int& iy, const int& iz, const int& it) const;
    void*                pixels(void);
    int                  type(void) const;

private:
    // Private methods
    void  init_members(void);
    void  copy_members(const GFitsImageByte& image);
    void  free_members(void);
    void  alloc_data(void);
    void  init_data(void);
    void  release_data(void);
    void  construct_data(const unsigned char* pixels);
    void  load_data(void) const;
    void  alloc_nulval(const void* value);
    void* ptr_data(void) { return m_pixels; }
    void* ptr_nulval(void) { return m_nulval; }

    // Private data area
    unsigned char* m_pixels;      //!< Pixels
    unsigned char* m_nulval;      //!< NULL value
};

#endif /* GFITSIMAGEBYTE_HPP */
