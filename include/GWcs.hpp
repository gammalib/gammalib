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
    GWcs(const GWcs& wcs);
    GWcs(GSkyDir& crval, const double& crpix1, const double& crpix2,
         const double& cdelt1, const double& cdelt2,
         const std::string& coords);
    virtual ~GWcs(void);

    // Operators
    virtual GWcs& operator= (const GWcs& wcs);

    // Virtual Methods
    virtual std::string type(void) const = 0;
    virtual void        read(const GFitsHDU* hdu) = 0;
    virtual void        write(GFitsHDU* hdu) = 0;
    virtual GSkyDir     pix2dir(const int& pix) = 0;
    virtual int         dir2pix(GSkyDir dir) const = 0;
    virtual double      omega(const int& pix) const = 0;

    // Implemented methods
    std::string coordsys(void) const;
    void        coordsys(const std::string& coordsys);

private:
    // Private methods
    void          init_members(void);
    void          copy_members(const GWcs& wcs);
    void          free_members(void);
    virtual GWcs* clone(void) const = 0;

protected:
    // Protected data area
    int          m_coordsys;     //!< 0=celestial, 1=galactic

    // WCS projection parameters
    std::string  m_ctype[2];     //!< Coordinate strings
    double       m_cd[4];        //!< Astrometry parameters
    double       m_cdelt[2];     //!< Increment at reference point (deg/pixel)
    double       m_crpix[2];     //!< x and y coordinates of the reference point
    double       m_crval[2];     //!< Coordinates of reference point
    double       m_longpole;     //!< Native longitude of North Pole
    double       m_latpole;      //!< Native latitude of North Pole
    double       m_pv2[2];       //!< Additional parameters
};

#endif /* GWCS_HPP */
