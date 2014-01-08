/***************************************************************************
 *          GSkyProjection.hpp - Abstract sky projection base class        *
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
 * @file GSkyProjection.hpp
 * @brief Abstract sky projection base class definition
 * @author Juergen Knoedlseder
 */

#ifndef GSKYPROJECTION_HPP
#define GSKYPROJECTION_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GFitsHDU.hpp"
#include "GSkyDir.hpp"
#include "GSkyPixel.hpp"


/***********************************************************************//**
 * @class GSkyProjection
 *
 * @brief Abstract sky projection base class
 *
 * This class defines an abstract projection from sky coordinates into
 * pixel coordinates. Sky coordinates are implemented using the GSkyDir
 * class, pixel coordinates are implemented using the GSkyPixel class.
 ***************************************************************************/
class GSkyProjection : public GBase {

    // Operator friends
    friend bool operator==(const GSkyProjection &a, const GSkyProjection &b);
    friend bool operator!=(const GSkyProjection &a, const GSkyProjection &b);

public:
    // Constructors and destructors
    GSkyProjection(void);
    GSkyProjection(const GSkyProjection& proj);
    virtual ~GSkyProjection(void);

    // Operators
    virtual GSkyProjection& operator=(const GSkyProjection& proj);

    // Pure virtual methods
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
    virtual std::string     print(const GChatter& chatter = NORMAL) const = 0;

    // Virtual methods
    virtual std::string coordsys(void) const;
    virtual void        coordsys(const std::string& coordsys);

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GSkyProjection& proj);
    void         free_members(void);
    virtual bool compare(const GSkyProjection& proj) const = 0;

    // Protected members
    int m_coordsys;   //!< 0=EQU, 1=GAL
};

#endif /* GSKYPROJECTION_HPP */
