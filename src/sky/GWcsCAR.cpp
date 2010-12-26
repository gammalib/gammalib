/***************************************************************************
 *                 GWcsCAR.cpp  -  Cartesian projection class              *
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
#define G_XY2DIR                                      "GWcsCAR::pix2dir(int)"
#define G_DIR2XY                                  "GWcsCAR::dir2pix(GSkyDir)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Local prototypes ___________________________________________________ */

/* __ Constants __________________________________________________________ */

/* __ Static conversion arrays ___________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
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

    // Initialise derived projection parameters
    wcs_init(0.0);  // theta0 = 0.0

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor from FITS HDU table
 *
 * @param[in] hdu FITS HDU.
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
 * @param[in] wcs World Coordinate System.
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
GWcsCAR::~GWcsCAR(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] wcs World Coordinate System.
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
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GWcsCAR::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GWcs::free_members();

    // Initialise members
    this->GWcs::init_members();
    init_members();

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


/***********************************************************************//**
 * @brief Read WCS definiton from FITS header
 *
 * @param[in] hdu FITS HDU.
 ***************************************************************************/
void GWcsCAR::read(const GFitsHDU* hdu)
{
    // Clear object
    clear();

    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Read standard WCS definition
        wcs_read(hdu);

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write Healpix definiton into FITS HDU
 *
 * @param[in] hdu FITS HDU.
 ***************************************************************************/
void GWcsCAR::write(GFitsHDU* hdu) const
{
    // Continue only if HDU pointer is valid
    if (hdu != NULL) {

        // Write standard WCS definition
        wcs_write(hdu);

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns solid angle of pixel in units of steradians
 *
 * @param[in] pix Pixel index (x,y)
 *
 * Estimate solid angles of pixel by compuing the coordinates in the 4 pixel
 * corners. The surface is computed using a cartesian approximation:
 *           a
 *     1-----------2                 a+b
 * h  /             \    where A = h ---
 *   4---------------3                2
 *           b
 * This is a brute force technique that works sufficiently well for non-
 * rotated sky maps. Something more intelligent should be implemented in
 * the future.
 *
 * @todo Implement accurate solid angle computation (so far only brute force
 *       estimation)
 ***************************************************************************/
double GWcsCAR::omega(const GSkyPixel& pix) const
{
    // Bypass const correctness
    GWcsCAR* ptr = (GWcsCAR*)this;

    // Get the sky directions of the 4 corners
    GSkyDir dir1 = ptr->xy2dir(GSkyPixel(pix.x()-0.5, pix.y()-0.5));
    GSkyDir dir2 = ptr->xy2dir(GSkyPixel(pix.x()+0.5, pix.y()-0.5));
    GSkyDir dir3 = ptr->xy2dir(GSkyPixel(pix.x()+0.5, pix.y()+0.5));
    GSkyDir dir4 = ptr->xy2dir(GSkyPixel(pix.x()-0.5, pix.y()+0.5));
    GSkyDir dir5 = ptr->xy2dir(GSkyPixel(pix.x(), pix.y()-0.5));
    GSkyDir dir6 = ptr->xy2dir(GSkyPixel(pix.x(), pix.y()+0.5));

    // Compute distances between sky directions
    double a = dir1.dist(dir2);
    double b = dir3.dist(dir4);
    double h = dir5.dist(dir6);

    // Compute solid angle
    double omega = 0.5*(h*(a+b));

    // Return solid angle
    return omega;
}


/***********************************************************************//**
 * @brief Print WCS information
 ***************************************************************************/
std::string GWcsCAR::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GWcsCAR ===\n");
    result.append(wcs_dump());

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GWcsCAR::init_members(void)
{
    // Initialise members
    m_type = "CAR";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] wcs World Coordinate System.
 ***************************************************************************/
void GWcsCAR::copy_members(const GWcsCAR& wcs)
{
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
 * @brief Projects from standard to native coordinates
 *
 * @param[in,out] coord Pointer to coordinate vector.
 *
 * Transforms coordinate vector from the standard coordinate system to the
 * native coordinate system.
 ***************************************************************************/
void GWcsCAR::std2nat(GVector *coord) const
{
    // Perform CAR projection
    *coord *= rad2deg;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Projects from native to standard coordinates
 *
 * @param[in,out] coord Pointer to coordinate vector.
 *
 * Transforms coordinate vector from the native coordinate system to the
 * standard coordinate system.
 ***************************************************************************/
void GWcsCAR::nat2std(GVector *coord) const
{
    // Perform CAR projection
    *coord *= deg2rad;

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
