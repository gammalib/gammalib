/***************************************************************************
 *                GFitsImage.i - Abstract FITS image base class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2020 by Juergen Knoedlseder                         *
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
 * @file GFitsImage.i
 * @brief Abstract FITS image base class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsImage.hpp"
#include "GFitsImageByte.hpp"
#include "GFitsImageDouble.hpp"
#include "GFitsImageFloat.hpp"
#include "GFitsImageLong.hpp"
#include "GFitsImageLongLong.hpp"
#include "GFitsImageSByte.hpp"
#include "GFitsImageShort.hpp"
#include "GFitsImageULong.hpp"
#include "GFitsImageUShort.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GFitsImage
 *
 * @brief Abstract FITS image base class
 ***************************************************************************/
class GFitsImage : public GFitsHDU {
public:
    // Constructors and destructors
    GFitsImage(void);
    GFitsImage(const int& bitpix, const int& nx);
    GFitsImage(const int& bitpix, const int& nx, const int& ny);
    GFitsImage(const int& bitpix, const int& nx, const int& ny, const int& nz);
    GFitsImage(const int& bitpix, const int& nx, const int& ny, const int& nz, const int& nt);
    GFitsImage(const int& bitpix, const std::vector<int>& naxes);
    GFitsImage(const GFitsImage& image);
    virtual ~GFitsImage(void);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GFitsImage* clone(void) const = 0;
    virtual std::string classname(void) const = 0;
    virtual double      pixel(const int& ix) const = 0;
    virtual double      pixel(const int& ix, const int& iy) const = 0;
    virtual double      pixel(const int& ix, const int& iy, const int& iz) const = 0;
    virtual double      pixel(const int& ix, const int& iy, const int& iz, const int& it) const = 0;
    virtual int         type(void) const = 0;

    // Implemented pure virtual methods
    HDUType exttype(void) const;

    // Base class methods
    const int&  npix(void) const;
    const int&  bitpix(void) const;
    const int&  naxis(void) const;
    int         naxes(const int& axis) const;
    const int&  anynul(void) const;
    void        nulval(const void* value);
    const void* nulval(void) const;
};


/***********************************************************************//**
 * @brief GFitsImage class extension
 ***************************************************************************/
%extend GFitsImage {
};
