/***************************************************************************
 *       GFitsImageLongLong.hpp - Long long integer FITS image class       *
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
 * @file GFitsImageLongLong.hpp
 * @brief Long long integer FITS image class definition
 * @author Juergen Knoedlseder
 */

#ifndef GFITSIMAGELONGLONG_HPP
#define GFITSIMAGELONGLONG_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsImage.hpp"


/***********************************************************************//**
 * @class GFitsImageLongLong
 *
 * @brief Long long integer FITS image class
 ***************************************************************************/
class GFitsImageLongLong : public GFitsImage {

public:
    // Constructors and destructors
    GFitsImageLongLong(void);
    GFitsImageLongLong(const int& nx, const long long* pixels = NULL);
    GFitsImageLongLong(const int& nx, const int& ny, const long long* pixels = NULL);
    GFitsImageLongLong(const int& nx, const int& ny, const int& nz, const long long* pixels = NULL);
    GFitsImageLongLong(const int& nx, const int& ny, const int& nz, const int& nt, const long long* pixels = NULL);
    GFitsImageLongLong(const std::vector<int>& naxes, const long long* pixels = NULL);
    GFitsImageLongLong(const GFitsImageLongLong& image);
    virtual ~GFitsImageLongLong(void);

    // Operators
    GFitsImageLongLong& operator=(const GFitsImageLongLong& image);
    long long&          operator()(const int& ix);
    long long&          operator()(const int& ix, const int& iy);
    long long&          operator()(const int& ix, const int& iy, const int& iz);
    long long&          operator()(const int& ix, const int& iy, const int& iz, const int& it);
    const long long&    operator()(const int& ix) const;
    const long long&    operator()(const int& ix, const int& iy) const;
    const long long&    operator()(const int& ix, const int& iy, const int& iz) const;
    const long long&    operator()(const int& ix, const int& iy, const int& iz, const int& it) const;

    // Methods
    void                clear(void);
    GFitsImageLongLong* clone(void) const;
    std::string         classname(void) const;
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
    int                 type(void) const;

private:
    // Private methods
    void  init_members(void);
    void  copy_members(const GFitsImageLongLong& image);
    void  free_members(void);
    void  alloc_data(void);
    void  init_data(void);
    void  release_data(void);
    void  construct_data(const long long* pixels);
    void  load_data(void) const;
    void  alloc_nulval(const void* value);
    void* ptr_data(void) { return m_pixels; }
    void* ptr_nulval(void) { return m_nulval; }

    // Private data area
    long long* m_pixels;      //!< Pixels
    long long* m_nulval;      //!< NULL value
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GFitsImageLongLong").
 ***************************************************************************/
inline
std::string GFitsImageLongLong::classname(void) const
{
    return ("GFitsImageLongLong");
}

#endif /* GFITSIMAGELONGLONG_HPP */
