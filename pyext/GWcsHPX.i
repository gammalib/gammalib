/***************************************************************************
 *                  GWcsHPX.i  -  Healpix projection class                 *
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
 * @file GWcsHPX.i
 * @brief HealPix projection class Python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GWcsHPX.hpp"
%}


/***********************************************************************//**
 * @class GWcsHPX
 *
 * @brief HealPix projection class Python interface defintion
 ***************************************************************************/
class GWcsHPX : public GWcs {
public:
    // Constructors and destructors
    GWcsHPX(void);
    explicit GWcsHPX(const int& nside, const std::string& ordering = "NESTED",
                     const std::string& coordsys = "GAL");
    explicit GWcsHPX(const GFitsHDU* hdu);
    GWcsHPX(const GWcsHPX& wcs);
    virtual ~GWcsHPX(void);

    // Implemented pure virtual methods
    virtual void        clear(void);
    virtual GWcsHPX*    clone(void) const;
    virtual std::string code(void) const;
    virtual std::string name(void) const;
    virtual void        read(const GFitsHDU* hdu);
    virtual void        write(GFitsHDU* hdu) const;
    virtual double      omega(const int& pix) const;
    virtual double      omega(const GSkyPixel& pix) const;
    virtual GSkyDir     pix2dir(const int& pix) const;
    virtual int         dir2pix(const GSkyDir& dir) const;
    virtual GSkyDir     xy2dir(const GSkyPixel& pix) const;
    virtual GSkyPixel   dir2xy(const GSkyDir& dir) const;

    // Additional class specific methods
    int          npix(void) const;
    int          nside(void) const;
    std::string  ordering(void) const;
    void         ordering(const std::string& ordering);
};


/***********************************************************************//**
 * @brief GWcsHPX class extension
 ***************************************************************************/
%extend GWcsHPX {
    GWcsHPX copy() {
        return (*self);
    }
};
