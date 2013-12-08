/***************************************************************************
 *                   GHealpix.i - Healpix projection class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @file GHealpix.i
 * @brief HealPix projection class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GHealpix.hpp"
%}


/***********************************************************************//**
 * @class GHealpix
 *
 * @brief HealPix projection class
 ***************************************************************************/
class GHealpix : public GSkyProjection {
public:
    // Constructors and destructors
    GHealpix(void);
    explicit GHealpix(const int& nside, const std::string& ordering = "NESTED",
                      const std::string& coordsys = "GAL");
    explicit GHealpix(const GFitsHDU& hdu);
    GHealpix(const GHealpix& wcs);
    virtual ~GHealpix(void);

    // Implemented pure virtual methods
    virtual void        clear(void);
    virtual GHealpix*   clone(void) const;
    virtual int         size(void) const;
    virtual std::string code(void) const;
    virtual std::string name(void) const;
    virtual void        read(const GFitsHDU& hdu);
    virtual void        write(GFitsHDU& hdu) const;
    virtual double      solidangle(const GSkyPixel& pixel) const;
    virtual GSkyDir     pix2dir(const GSkyPixel& pixel) const;
    virtual GSkyPixel   dir2pix(const GSkyDir& dir) const;

    // Other methods
    const int&   npix(void) const;
    const int&   nside(void) const;
    std::string  ordering(void) const;
    void         ordering(const std::string& ordering);
};


/***********************************************************************//**
 * @brief GHealpix class extension
 ***************************************************************************/
%extend GHealpix {
    GHealpix copy() {
        return (*self);
    }
};
