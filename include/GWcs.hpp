/***************************************************************************
 *          GWcs.hpp  -  World Coordinate System virtual base class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
 * @file GWcs.hpp
 * @brief World Coordinate System virtual base class definition
 * @author Juergen Knoedlseder
 */

#ifndef GWCS_HPP
#define GWCS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GBase.hpp"
#include "GLog.hpp"
#include "GFitsHDU.hpp"
#include "GSkyDir.hpp"
#include "GSkyPixel.hpp"


/***********************************************************************//**
 * @class GWcs
 *
 * @brief GWcs virtual base class interface defintion
 ***************************************************************************/
class GWcs : public GBase {

    // Operator friends
    friend bool operator== (const GWcs &a, const GWcs &b);
    friend bool operator!= (const GWcs &a, const GWcs &b);

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GWcs& wcs);
    friend GLog&         operator<< (GLog& log, const GWcs& wcs);

public:
    // Constructors and destructors
    GWcs(void);
    GWcs(const GWcs& wcs);
    virtual ~GWcs(void);

    // Operators
    virtual GWcs& operator= (const GWcs& wcs);

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
    virtual std::string print(void) const = 0;

    // Virtual methods
    virtual std::string coordsys(void) const;
    virtual void        coordsys(const std::string& coordsys);

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GWcs& wcs);
    void         free_members(void);
    virtual bool compare(const GWcs& wcs) const = 0;

    // Protected members
    int m_coordsys;   //!< 0=celestial, 1=galactic, 2=ecl, 3=hel, 4=sgl
};

#endif /* GWCS_HPP */
