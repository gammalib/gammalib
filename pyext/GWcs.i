/***************************************************************************
 *          GWcs.i  -  World Coordinate System virtual base class          *
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
 * @file GWcs.i
 * @brief World Coordinate System virtual base class Python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GWcs.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GWcs
 *
 * @brief World Coordinate System virtual base class Python interface
 ***************************************************************************/
class GWcs {
public:
    // Constructors and destructors
    GWcs(void);
    GWcs(const GWcs& wcs);
    virtual ~GWcs(void);

    // Pure virtual methods (not implemented)
    virtual void        clear(void) = 0;
    virtual GWcs*       clone(void) const = 0;
    virtual std::string code(void) const = 0;
    virtual std::string name(void) const = 0;
    virtual void        read(const GFitsHDU* hdu) = 0;
    virtual void        write(GFitsHDU* hdu) const = 0;
    virtual double      omega(const int& pix) const = 0;
    virtual double      omega(const GSkyPixel& pix) const = 0;
    virtual GSkyDir     pix2dir(const int& pix) const = 0;
    virtual int         dir2pix(const GSkyDir& dir) const = 0;
    virtual GSkyDir     xy2dir(const GSkyPixel& pix) const = 0;
    virtual GSkyPixel   dir2xy(const GSkyDir& dir) const = 0;

    // Virtual methods
    virtual std::string coordsys(void) const;
    virtual void        coordsys(const std::string& coordsys);
};




/***********************************************************************//**
 * @brief GWcs class extension
 ***************************************************************************/
%extend GWcs {
    char *__str__() {
        return tochar(self->print());
    }
    bool __is__(const GWcs &a) {
            return (*self) == a;
    }
};
