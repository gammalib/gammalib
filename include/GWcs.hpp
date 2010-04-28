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
    // Protected data area
    int     m_coordsys;     //!< 0=celestial, 1=galactic
    double  m_crval[2];     //!< Coordinates of reference point
    double  m_crpix[2];     //!< x and y coordinates of the reference point
    double  m_cdelt[2];     //!< Increment at reference point (deg/pixel)
    double  m_longpole;     //!< Native longitude of North Pole
    double  m_latpole;      //!< Native latitude of North Pole
    GMatrix m_cd;           //!< Astrometry parameters (2x2 matrix, deg/pixel)
    GMatrix m_invcd;        //!< Inverse of CD matrix
    GVector m_pv2;          //!< Projection parameters (up to 21)
};

#endif /* GWCS_HPP */
