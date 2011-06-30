/***************************************************************************
 *          GWcs.hpp  -  World Coordinate System virtual base class        *
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
 * @file GWcs.hpp
 * @brief GWcs virtual base class definition.
 * @author J. Knodlseder
 */

#ifndef GWCS_HPP
#define GWCS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GLog.hpp"
#include "GFitsHDU.hpp"
#include "GSkyDir.hpp"
#include "GSkyPixel.hpp"
#include "GMatrix.hpp"
#include "GVector.hpp"


/***********************************************************************//**
 * @class GWcs
 *
 * @brief GWcs virtual base class interface defintion
 ***************************************************************************/
class GWcs {

    // Friend classes
    friend class GSkymap;

    // Operator friends
    friend bool operator== (const GWcs &a, const GWcs &b);
    friend bool operator!= (const GWcs &a, const GWcs &b);

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GWcs& wcs);
    friend GLog&         operator<< (GLog& log, const GWcs& wcs);

public:
    // Constructors and destructors
    GWcs(void);
    //explicit GWcs(const std::string& coords,
    //              const double& crval1, const double& crval2,
    //              const double& crpix1, const double& crpix2,
    //              const double& cdelt1, const double& cdelt2);
    GWcs(const GWcs& wcs);
    virtual ~GWcs(void);

    // Operators
    virtual GWcs& operator= (const GWcs& wcs);

    // Pure virtual methods (not implemented)
    virtual void        clear(void) = 0;
    virtual GWcs*       clone(void) const = 0;
    virtual void        read(const GFitsHDU* hdu) = 0;
    virtual void        write(GFitsHDU* hdu) const = 0;
    virtual std::string print(void) const = 0;

    // Virtual methods
    virtual std::string type(void) const;
    virtual std::string coordsys(void) const;
    virtual void        coordsys(const std::string& coordsys);
    virtual double      omega(const int& pix) const;
    virtual double      omega(const GSkyPixel& pix) const;
    virtual GSkyDir     pix2dir(const int& pix) const;
    virtual int         dir2pix(GSkyDir dir) const;
    virtual GSkyDir     xy2dir(const GSkyPixel& pix) const;
    virtual GSkyPixel   dir2xy(GSkyDir dir) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GWcs& wcs);
    void free_members(void);

    // Protected methods (providing services to derived classes)
    virtual void std2nat(GVector *coord) const = 0;
    virtual void nat2std(GVector *coord) const = 0;
    GVector      wcs_dir2native(GSkyDir dir) const;
    GSkyDir      wcs_native2dir(GVector native) const;
    void         wcs_init(const double& theta0);
    GVector      wcs_getpole(const double& theta0);
    GMatrix      wcs_get_rot(void);
    //void         wcs_read(const GFitsHDU* hdu);
    //void         wcs_write(GFitsHDU* hdu) const;
    std::string  wcs_crval1(void) const;
    std::string  wcs_crval2(void) const;
    //std::string  wcs_dump(void) const;

    // Protected members (astr structure)
    std::string m_type;     //!< WCS type
    int         m_coordsys; //!< 0=celestial, 1=galactic
    int         m_reverse;  //!< Reverse axes (1=true)
    GVector     m_crval;    //!< Value of reference point
    GVector     m_crpix;    //!< Pixel of reference point (starting from 1)
    GVector     m_cdelt;    //!< Increment at reference point (deg/pixel)
    GVector     m_npole;    //!< Native coordinate of North Pole
    GMatrix     m_cd;       //!< Astrometry parameters (2x2 matrix, deg/pixel)
    GVector     m_pv2;      //!< Projection parameters (up to 21)

    // Derived members
    //double  m_theta0;       //!< Native latitude of the fiducial point
    GVector m_refval;       //!< Ordered value of reference point
    GVector m_refpix;       //!< Pixel of reference point (starting from 0)
    GMatrix m_invcd;        //!< Inverse of CD matrix
    GVector m_native_pole;  //!< Coordinates of native pole
    GMatrix m_rot;          //!< Rotation matrix
    GMatrix m_trot;         //!< Transpose of rotation matrix
};

#endif /* GWCS_HPP */
