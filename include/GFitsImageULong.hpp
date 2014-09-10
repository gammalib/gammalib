/***************************************************************************
 *          GFitsImageULong.hpp - Unsigned long image FITS class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
 * @file GFitsImageULong.hpp
 * @brief Unsigned long image FITS class definition
 * @author Juergen Knoedlseder
 */

#ifndef GFITSIMAGEULONG_HPP
#define GFITSIMAGEULONG_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsImage.hpp"


/***********************************************************************//**
 * @class GFitsImageULong
 *
 * @brief Unsigned long image FITS class
 ***************************************************************************/
class GFitsImageULong : public GFitsImage {

public:
    // Constructors and destructors
    GFitsImageULong(void);
    GFitsImageULong(const int& nx, const unsigned long* pixels = NULL);
    GFitsImageULong(const int& nx, const int& ny, const unsigned long* pixels = NULL);
    GFitsImageULong(const int& nx, const int& ny, const int& nz, const unsigned long* pixels = NULL);
    GFitsImageULong(const int& nx, const int& ny, const int& nz, const int& nt, const unsigned long* pixels = NULL);
    GFitsImageULong(const int& naxis, const int* naxes, const unsigned long* pixels = NULL);
    GFitsImageULong(const GFitsImageULong& image);
    virtual ~GFitsImageULong(void);

    // Operators
    GFitsImageULong&     operator=(const GFitsImageULong& image);
    unsigned long&       operator()(const int& ix);
    unsigned long&       operator()(const int& ix, const int& iy);
    unsigned long&       operator()(const int& ix, const int& iy, const int& iz);
    unsigned long&       operator()(const int& ix, const int& iy, const int& iz, const int& it);
    const unsigned long& operator()(const int& ix) const;
    const unsigned long& operator()(const int& ix, const int& iy) const;
    const unsigned long& operator()(const int& ix, const int& iy, const int& iz) const;
    const unsigned long& operator()(const int& ix, const int& iy, const int& iz, const int& it) const;

    // Methods
    void                 clear(void);
    GFitsImageULong*     clone(void) const;
    std::string          classname(void) const;
    unsigned long&       at(const int& ix);
    unsigned long&       at(const int& ix, const int& iy);
    unsigned long&       at(const int& ix, const int& iy, const int& iz);
    unsigned long&       at(const int& ix, const int& iy, const int& iz, const int& it);
    const unsigned long& at(const int& ix) const;
    const unsigned long& at(const int& ix, const int& iy) const;
    const unsigned long& at(const int& ix, const int& iy, const int& iz) const;
    const unsigned long& at(const int& ix, const int& iy, const int& iz, const int& it) const;
    double               pixel(const int& ix) const;
    double               pixel(const int& ix, const int& iy) const;
    double               pixel(const int& ix, const int& iy, const int& iz) const;
    double               pixel(const int& ix, const int& iy, const int& iz, const int& it) const;
    void*                pixels(void);
    int                  type(void) const;

private:
    // Private methods
    void  init_members(void);
    void  copy_members(const GFitsImageULong& image);
    void  free_members(void);
    void  alloc_data(void);
    void  init_data(void);
    void  release_data(void);
    void  construct_data(const unsigned long* pixels);
    void  load_data(void) const;
    void  alloc_nulval(const void* value);
    void* ptr_data(void) { return m_pixels; }
    void* ptr_nulval(void) { return m_nulval; }

    // Private data area
    unsigned long* m_pixels;      //!< Pixels
    unsigned long* m_nulval;      //!< NULL value
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GFitsImageULong").
 ***************************************************************************/
inline
std::string GFitsImageULong::classname(void) const
{
    return ("GFitsImageULong");
}

#endif /* GFITSIMAGEULONG_HPP */
