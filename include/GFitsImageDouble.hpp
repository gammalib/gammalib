/***************************************************************************
 *        GFitsImageDouble.hpp - Double precision FITS image class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2016 by Juergen Knoedlseder                         *
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
 * @brief Double precision FITS image class definition
 * @author Juergen Knoedlseder
 */

#ifndef GFITSIMAGEDOUBLE_HPP
#define GFITSIMAGEDOUBLE_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsImage.hpp"


/***********************************************************************//**
 * @class GFitsImageDouble
 *
 * @brief Double precision FITS image class
 ***************************************************************************/
class GFitsImageDouble : public GFitsImage {

public:
    // Constructors and destructors
    GFitsImageDouble(void);
    GFitsImageDouble(const int& nx, const double* pixels = NULL);
    GFitsImageDouble(const int& nx, const int& ny, const double* pixels = NULL);
    GFitsImageDouble(const int& nx, const int& ny, const int& nz, const double* pixels = NULL);
    GFitsImageDouble(const int& nx, const int& ny, const int& nz, const int& nt, const double* pixels = NULL);
    GFitsImageDouble(const std::vector<int>& naxes, const double* pixels = NULL);
    GFitsImageDouble(const GFitsImageDouble& image);
    virtual ~GFitsImageDouble(void);

    // Operators
    GFitsImageDouble& operator=(const GFitsImageDouble& image);
    double&           operator()(const int& ix);
    double&           operator()(const int& ix, const int& iy);
    double&           operator()(const int& ix, const int& iy, const int& iz);
    double&           operator()(const int& ix, const int& iy, const int& iz, const int& it);
    const double&     operator()(const int& ix) const;
    const double&     operator()(const int& ix, const int& iy) const;
    const double&     operator()(const int& ix, const int& iy, const int& iz) const;
    const double&     operator()(const int& ix, const int& iy, const int& iz, const int& it) const;

    // Methods
    void              clear(void);
    GFitsImageDouble* clone(void) const;
    std::string       classname(void) const;
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


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GFitsImageDouble").
 ***************************************************************************/
inline
std::string GFitsImageDouble::classname(void) const
{
    return ("GFitsImageDouble");
}

#endif /* GFITSIMAGEDOUBLE_HPP */
