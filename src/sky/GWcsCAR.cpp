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

    // Initialise derived projection parameters
    wcs_init(0.0);  // theta0 = 0.0

    // Set projection function pointers
    m_std2nat = (_wcspf)&GWcsCAR::std2nat;
    m_nat2std = (_wcspf)&GWcsCAR::nat2std;

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
 * @brief Read WCS definiton from FITS header
 *
 * @param[in] hdu FITS HDU containing the Healpix definition.
 ***************************************************************************/
void GWcsCAR::read(const GFitsHDU* hdu)
{
    // Free memory and initialise members
    free_members();
    init_members();

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
 * @param[in] hdu FITS HDU to which the Healpix definition will be written.
 ***************************************************************************/
void GWcsCAR::write(GFitsHDU* hdu) const
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Write standard WCS definition
        wcs_write(hdu);

    } // endif: HDU was valid

    // Return
    return;
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
    m_type = "CAR";

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
    wcs.wcs_dump(os);

    // Return output stream
    return os;
}
