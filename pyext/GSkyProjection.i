/***************************************************************************
 *           GSkyProjection.i - Abstract sky projection base class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GSkyProjection.i
 * @brief Abstract sky projection base class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSkyProjection.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GSkyProjection
 *
 * @brief Abstract sky projection base class
 *
 * This class defines an abstract projection from sky coordinates into
 * pixel coordinates.
 ***************************************************************************/
class GSkyProjection : public GBase {
public:
    // Constructors and destructors
    GSkyProjection(void);
    GSkyProjection(const GSkyProjection& proj);
    virtual ~GSkyProjection(void);

    // Pure virtual methods (not implemented)
    virtual void            clear(void) = 0;
    virtual GSkyProjection* clone(void) const = 0;
    virtual int             size(void) const = 0;
    virtual std::string     code(void) const = 0;
    virtual std::string     name(void) const = 0;
    virtual void            read(const GFitsHDU& hdu) = 0;
    virtual void            write(GFitsHDU& hdu) const = 0;
    virtual double          solidangle(const GSkyPixel& pixel) const = 0;
    virtual GSkyDir         pix2dir(const GSkyPixel& pixel) const = 0;
    virtual GSkyPixel       dir2pix(const GSkyDir& dir) const = 0;

    // Virtual methods
    virtual std::string coordsys(void) const;
    virtual void        coordsys(const std::string& coordsys);
};


/***********************************************************************//**
 * @brief GSkyProjection class extension
 ***************************************************************************/
%extend GSkyProjection {
    bool __is__(const GSkyProjection &proj) {
            return (*self) == proj;
    }
};
