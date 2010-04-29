/***************************************************************************
 *          GWcs.hpp  -  World Coordinate System virtual base class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
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

public:
    // Constructors and destructors
    GWcs(void);
    GWcs(const std::string& coords,
         const double& crval1, const double& crval2,
         const double& crpix1, const double& crpix2,
         const double& cdelt1, const double& cdelt2);
    GWcs(const GWcs& wcs);
    virtual ~GWcs(void);

    // Operators
    virtual GWcs& operator= (const GWcs& wcs);

    // Pure virtual methods (not implemented)
    virtual std::string type(void) const = 0;
    virtual GSkyDir     pix2dir(const int& pix) = 0;
    virtual GSkyDir     xy2dir(const GSkyPixel& pix) = 0;
    virtual int         dir2pix(GSkyDir dir) const = 0;
    virtual GSkyPixel   dir2xy(GSkyDir dir) const = 0;
    virtual double      omega(const int& pix) const = 0;
    virtual double      omega(const GSkyPixel& pix) const = 0;

    // Virtual methods (implemented)
    virtual void        read(const GFitsHDU* hdu);
    virtual void        write(GFitsHDU* hdu);
    virtual std::string coordsys(void) const;
    virtual void        coordsys(const std::string& coordsys);

private:
    // Private methods
    void          init_members(void);
    void          copy_members(const GWcs& wcs);
    void          free_members(void);
    virtual GWcs* clone(void) const = 0;

protected:
    // Protected methods
    GVector wcs_dir2native(GSkyDir dir) const;
    GSkyDir wcs_native2dir(GVector native) const;
    GVector wcs_getpole(const double& theta0);
    GMatrix wcs_get_rot(void);
    void    dump_wcs(std::ostream& os) const;
    
    // Protected data area (astr structure)
    int     m_coordsys;     //!< 0=celestial, 1=galactic
    int     m_reverse;      //!< Reverse axes (1=true)
    GVector m_crval;        //!< Value of reference point
    GVector m_crpix;        //!< Pixel of reference point (starting from 1)
    GVector m_cdelt;        //!< Increment at reference point (deg/pixel)
    GVector m_npole;        //!< Native coordinate of North Pole
    GMatrix m_cd;           //!< Astrometry parameters (2x2 matrix, deg/pixel)
    GVector m_pv2;          //!< Projection parameters (up to 21)
    
    // Derived parameters
    double  m_theta0;       //!< Native latitude of the fiducial point
    GVector m_refval;       //!< Ordered value of reference point
    GVector m_refpix;       //!< Pixel of reference point (starting from 0)
    GMatrix m_invcd;        //!< Inverse of CD matrix
    GVector m_native_pole;  //!< Coordinates of native pole
    GMatrix m_rot;          //!< Rotation matrix
    GMatrix m_trot;         //!< Transpose of rotation matrix
};

#endif /* GWCS_HPP */
