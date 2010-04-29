/***************************************************************************
 *                 GWcsCAR.cpp  -  Healpix projection class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
/**
 * @file GWcsCAR.cpp
 * @brief GWcsCAR class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GWcsCAR.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_XY2DIR                                 "GWcsCAR::xy2dir(GSkyPixel)"
#define G_DIR2XY                                   "GWcsCAR::dir2xy(GSkyDir)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DIR2XY_DEBUG                                      // Debug dir2xy
//#define G_XY2DIR_DEBUG                                      // Debug xy2dir

/* __ Local prototypes ___________________________________________________ */

/* __ Constants __________________________________________________________ */

/* __ Static conversion arrays ___________________________________________ */


/*==========================================================================
 =                                                                         =
 =                      GWcsCAR constructors/destructors                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GWcsCAR::GWcsCAR(void) : GWcs()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] coords Coordinate system.
 * @param[in] crval1 X value of reference pixel.
 * @param[in] crval1 Y value of reference pixel.
 * @param[in] crpix1 X index of reference pixel.
 * @param[in] crpix2 Y index of reference pixel.
 * @param[in] cdelt1 Increment in x direction at reference pixel (deg).
 * @param[in] cdelt2 Increment in y direction at reference pixel (deg).
 * @param[in] cd Astrometry parameters (2x2 matrix, deg/pixel).
 * @param[in] pv2 Projection parameters (length WCS type dependent).
 ***************************************************************************/
GWcsCAR::GWcsCAR(const std::string& coords,
                 const double& crval1, const double& crval2,
                 const double& crpix1, const double& crpix2,
                 const double& cdelt1, const double& cdelt2,
                 const GMatrix& cd, const GVector& pv2) :
                 GWcs(coords, crval1, crval2, crpix1, crpix2, cdelt1, cdelt2)

{
    // Initialise class members
    init_members();

    // Derive projection parameters
    m_theta0      = 0.0;
    m_native_pole = wcs_getpole(m_theta0);
    m_rot         = wcs_get_rot();
    m_trot        = transpose(m_rot);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor from FITS HDU table
 *
 * @param[in] hdu Pointer to FITS HDU.
 ***************************************************************************/
GWcsCAR::GWcsCAR(const GFitsHDU* hdu) : GWcs()
{
    // Initialise class members
    init_members();

    // Read Healpix defintion from FITS HDU
    read(hdu);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] wcs GWcsCAR instance which should be used for construction.
 ***************************************************************************/
GWcsCAR::GWcsCAR(const GWcsCAR& wcs) : GWcs(wcs)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(wcs);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GWcsCAR::~GWcsCAR()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                            GWcsCAR operators                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] wcs GWcsHPX instance to be assigned.
 ***************************************************************************/
GWcsCAR& GWcsCAR::operator= (const GWcsCAR& wcs)
{
    // Execute only if object is not identical
    if (this != &wcs) {

        // Copy base class members
        this->GWcs::operator=(wcs);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(wcs);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                          GWcsCAR public methods                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Returns WCS type
 ***************************************************************************/
std::string GWcsCAR::type(void) const
{
    // Return Healix type
    return "CAR";
}


/***********************************************************************//**
 * @brief Returns sky direction of pixel (method should never be called)
 *
 * @param[in] pix Pixel number (0,1,...,m_num_pixels).
 ***************************************************************************/
GSkyDir GWcsCAR::pix2dir(const int& pix)
{
    // Throw error
    throw GException::wcs(G_XY2DIR, 
                          "WCS method not defined for CAR projection.");

    // Define constant direction
    const GSkyDir dir;

    // Return
    return dir;
}


/***********************************************************************//**
 * @brief Returns sky direction of pixel
 *
 * @param[in] pix Sky pixel.
 *
 * Note that pixel indices start from 0.
 ***************************************************************************/
GSkyDir GWcsCAR::xy2dir(const GSkyPixel& pix)
{
    // Set constants
    const int __permutation[2] = {1,0};

    // Determine offset
    GVector offset = GVector(pix.x(), pix.y());

    // Determine xy
    GVector xy = m_cd * (offset - m_refpix);

    // Swap result if coordinates are reversed
    if (m_reverse)
        xy = perm(xy, __permutation);

    // Perform CAR map projection
    GVector native = xy * deg2rad;

    // Get sky direction for native coordinates
    GSkyDir dir = wcs_native2dir(native);
    
    // Debug: Dump transformation steps
    #if defined(G_XY2DIR_DEBUG)
    std::cout << "xy2dir: pixel=" << offset << " xy=" << xy
              << " native=" << native << " dir=" << dir << std::endl;
    #endif

    // Return
    return dir;
}


/***********************************************************************//**
 * @brief Returns pixel for sky direction (method should never be called)
 *
 * @param[in] dir Sky direction.
 ***************************************************************************/
int GWcsCAR::dir2pix(GSkyDir dir) const
{
    // Throw error
    throw GException::wcs(G_DIR2XY, 
                          "WCS method not defined for CAR projection.");
    
    // Return
    return 0;
}


/***********************************************************************//**
 * @brief Returns pixel of sky direction
 *
 * @param[in] dir Sky direction.
 *
 * Note that pixel indices start from 0.
 ***************************************************************************/
GSkyPixel GWcsCAR::dir2xy(GSkyDir dir) const
{
    // Set constants
    const int __permutation[2] = {1,0};
        
    // Get native coordinates for sky direction in radians
    GVector native = wcs_dir2native(dir);
    
    // Perform CAR map projection
    GVector xy = native * rad2deg;
    
    // Swap result if coordinates are reversed
    if (m_reverse)
        xy = perm(xy, __permutation);
    
    // Determine pixel offset
    GVector offset = m_invcd * xy + m_refpix;

    // Set sky pixel
    GSkyPixel pixel(offset(0), offset(1));

    // Debug: Dump transformation steps
    #if defined(G_DIR2XY_DEBUG)
    std::cout << "dir2xy: dir=" << dir << " native=" << native
              << " xy=" << xy << " pixel=" << offset << std::endl;
    #endif

    // Return sky pixel
    return pixel;
}


/***********************************************************************//**
 * @brief Returns solid angle of pixel
 *
 * @param[in] pix Pixel number (0,1,...,m_num_pixels).
 ***************************************************************************/
double GWcsCAR::omega(const int& pix) const
{
    // TODO: Implement
    
    // Return solid angle
    return 0.0;
}


/***********************************************************************//**
 * @brief Returns solid angle of pixel
 *
 * @param[in] pix Pixel index (x,y)
 ***************************************************************************/
double GWcsCAR::omega(const GSkyPixel& pix) const
{
    // TODO: Implement
    
    // Return solid angle
    return 0.0;
}


/*==========================================================================
 =                                                                         =
 =                          GWcsCAR private methods                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GWcsCAR::init_members(void)
{
    // Initialise members

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] wcs GWcsCAR instance from which members should be copied.
 ***************************************************************************/
void GWcsCAR::copy_members(const GWcsCAR& wcs)
{
    // Copy attributes

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GWcsCAR::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 ***************************************************************************/
GWcsCAR* GWcsCAR::clone(void) const
{
    return new GWcsCAR(*this);
}


/*==========================================================================
 =                                                                         =
 =                              GWcsCAR friends                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] wcs Healpix WCS definition to put in output stream
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GWcsCAR& wcs)
{
    // Put header in stream
    os << "=== GWcsCAR ===" << std::endl;
    
    // Add WCS information
    wcs.dump_wcs(os);

    // Return output stream
    return os;
}
