/***************************************************************************
 *         GFitsImageUShort.hpp - Unsigned short FITS image class          *
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
 * @file GFitsImageUShort.hpp
 * @brief Unsigned short FITS image class definition
 * @author Juergen Knoedlseder
 */

#ifndef GFITSIMAGEUSHORT_HPP
#define GFITSIMAGEUSHORT_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsImage.hpp"


/***********************************************************************//**
 * @class GFitsImageUShort
 *
 * @brief Unsigned short FITS image class
 ***************************************************************************/
class GFitsImageUShort : public GFitsImage {

public:
    // Constructors and destructors
    GFitsImageUShort(void);
    GFitsImageUShort(const int& nx, const unsigned short* pixels = NULL);
    GFitsImageUShort(const int& nx, const int& ny, const unsigned short* pixels = NULL);
    GFitsImageUShort(const int& nx, const int& ny, const int& nz, const unsigned short* pixels = NULL);
    GFitsImageUShort(const int& nx, const int& ny, const int& nz, const int& nt, const unsigned short* pixels = NULL);
    GFitsImageUShort(const int& naxis, const int* naxes, const unsigned short* pixels = NULL);
    GFitsImageUShort(const GFitsImageUShort& image);
    virtual ~GFitsImageUShort(void);

    // Operators
    GFitsImageUShort&     operator=(const GFitsImageUShort& image);
    unsigned short&       operator()(const int& ix);
    unsigned short&       operator()(const int& ix, const int& iy);
    unsigned short&       operator()(const int& ix, const int& iy, const int& iz);
    unsigned short&       operator()(const int& ix, const int& iy, const int& iz, const int& it);
    const unsigned short& operator()(const int& ix) const;
    const unsigned short& operator()(const int& ix, const int& iy) const;
    const unsigned short& operator()(const int& ix, const int& iy, const int& iz) const;
    const unsigned short& operator()(const int& ix, const int& iy, const int& iz, const int& it) const;

    // Methods
    void                  clear(void);
    GFitsImageUShort*     clone(void) const;
    std::string           classname(void) const;
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
    int                   type(void) const;

private:
    // Private methods
    void  init_members(void);
    void  copy_members(const GFitsImageUShort& image);
    void  free_members(void);
    void  alloc_data(void);
    void  init_data(void);
    void  release_data(void);
    void  construct_data(const unsigned short* pixels);
    void  load_data(void) const;
    void  alloc_nulval(const void* value);
    void* ptr_data(void) { return m_pixels; }
    void* ptr_nulval(void) { return m_nulval; }

    // Private data area
    unsigned short* m_pixels;      //!< Pixels
    unsigned short* m_nulval;      //!< NULL value
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GFitsImageUShort").
 ***************************************************************************/
inline
std::string GFitsImageUShort::classname(void) const
{
    return ("GFitsImageUShort");
}

#endif /* GFITSIMAGEUSHORT_HPP */
