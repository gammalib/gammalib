/***************************************************************************
 *            GFitsImageLong.hpp - Long integer FITS image class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2016 by Juergen Knoedlseder                         *
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
 * @file GFitsImageLong.hpp
 * @brief Long integer FITS image class definition
 * @author Juergen Knoedlseder
 */

#ifndef GFITSIMAGELONG_HPP
#define GFITSIMAGELONG_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsImage.hpp"


/***********************************************************************//**
 * @class GFitsImageLong
 *
 * @brief Long integer FITS image class
 ***************************************************************************/
class GFitsImageLong : public GFitsImage {

public:
    // Constructors and destructors
    GFitsImageLong(void);
    GFitsImageLong(const int& nx, const long* pixels = NULL);
    GFitsImageLong(const int& nx, const int& ny, const long* pixels = NULL);
    GFitsImageLong(const int& nx, const int& ny, const int& nz, const long* pixels = NULL);
    GFitsImageLong(const int& nx, const int& ny, const int& nz, const int& nt, const long* pixels = NULL);
    GFitsImageLong(const std::vector<int>& naxes, const long* pixels = NULL);
    GFitsImageLong(const GFitsImageLong& image);
    virtual ~GFitsImageLong(void);

    // Operators
    GFitsImageLong& operator=(const GFitsImageLong& image);
    long&           operator()(const int& ix);
    long&           operator()(const int& ix, const int& iy);
    long&           operator()(const int& ix, const int& iy, const int& iz);
    long&           operator()(const int& ix, const int& iy, const int& iz, const int& it);
    const long&     operator()(const int& ix) const;
    const long&     operator()(const int& ix, const int& iy) const;
    const long&     operator()(const int& ix, const int& iy, const int& iz) const;
    const long&     operator()(const int& ix, const int& iy, const int& iz, const int& it) const;

    // Methods
    void            clear(void);
    GFitsImageLong* clone(void) const;
    std::string     classname(void) const;
    long&           at(const int& ix);
    long&           at(const int& ix, const int& iy);
    long&           at(const int& ix, const int& iy, const int& iz);
    long&           at(const int& ix, const int& iy, const int& iz, const int& it);
    const long&     at(const int& ix) const;
    const long&     at(const int& ix, const int& iy) const;
    const long&     at(const int& ix, const int& iy, const int& iz) const;
    const long&     at(const int& ix, const int& iy, const int& iz, const int& it) const;
    double          pixel(const int& ix) const;
    double          pixel(const int& ix, const int& iy) const;
    double          pixel(const int& ix, const int& iy, const int& iz) const;
    double          pixel(const int& ix, const int& iy, const int& iz, const int& it) const;
    void*           pixels(void);
    int             type(void) const;

private:
    // Private methods
    void  init_members(void);
    void  copy_members(const GFitsImageLong& image);
    void  free_members(void);
    void  alloc_data(void);
    void  init_data(void);
    void  release_data(void);
    void  construct_data(const long* pixels);
    void  load_data(void) const;
    void  alloc_nulval(const void* value);
    void* ptr_data(void) { return m_pixels; }
    void* ptr_nulval(void) { return m_nulval; }

    // Private data area
    long* m_pixels;      //!< Pixels
    long* m_nulval;      //!< NULL value
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GFitsImageLong").
 ***************************************************************************/
inline
std::string GFitsImageLong::classname(void) const
{
    return ("GFitsImageLong");
}

#endif /* GFITSIMAGELONG_HPP */
