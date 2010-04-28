/***************************************************************************
 *           GWcs.cpp  -  World Coordinate System virtual base class       *
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
 * @file GWcs.cpp
 * @brief GWcs virtual base class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GTools.hpp"
#include "GWcs.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_COORDSYS_SET                          "GWcs::coordsys(std::string)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Local prototypes ___________________________________________________ */

/* __ Constants __________________________________________________________ */

/* __ Static conversion arrays ___________________________________________ */

/*==========================================================================
 =                                                                         =
 =                       GWcs constructors/destructors                     =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GWcs::GWcs(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Standard 2D sky map constructor
 *
 * @param[in] coords Coordinate system.
 * @param[in] crval1 X value of reference pixel.
 * @param[in] crval1 Y value of reference pixel.
 * @param[in] crpix1 X index of reference pixel.
 * @param[in] crpix2 Y index of reference pixel.
 * @param[in] cdelt1 Increment in x direction at reference pixel (deg).
 * @param[in] cdelt2 Increment in y direction at reference pixel (deg).
 *
 * Construct standard 2D sky map from standard definition parameters. This
 * method
 ***************************************************************************/
GWcs::GWcs(const std::string& coords,
           const double& crval1, const double& crval2,
           const double& crpix1, const double& crpix2,
           const double& cdelt1, const double& cdelt2)

{
    // Initialise class members
    init_members();

    //TODO: Check parameters

    // Set coordinate system
    coordsys(coords);

    // Set parameters
    m_crval[0] = crval1;
    m_crval[1] = crval2;
    m_crpix[0] = crpix1;
    m_crpix[1] = crpix2;
    m_cdelt[0] = cdelt1;
    m_cdelt[1] = cdelt2;

    // Set standard CD without rotation is a unit matrix
    m_cd(0,0) = m_cd(1,1) = 1.0;
    m_cd(0,1) = m_cd(1,0) = 0.0;

    // Compute inverse CD matrix
    m_invcd = invert(m_cd); //NOTE: invert method not yet implemented

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param wcs GWcs instance which should be used for construction
 ***************************************************************************/
GWcs::GWcs(const GWcs& wcs)
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
GWcs::~GWcs(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                            GWcs operators                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] wcs GWcs instance to be assigned
 ***************************************************************************/
GWcs& GWcs::operator= (const GWcs& wcs)
{
    // Execute only if object is not identical
    if (this != &wcs) {

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
 =                           GWcs public methods                           =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Read WCS definiton from FITS header.
 *
 * @param[in] hdu FITS HDU containing the WCS definition.
 ***************************************************************************/
void GWcs::read(const GFitsHDU* hdu)
{
    // Free memory and initialise members
    free_members();
    init_members();

    //TODO: Implement WCS definition reading

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write WCS definiton into FITS HDU.
 *
 * @param[in] hdu FITS HDU into which WCS definition will be written.
 ***************************************************************************/
void GWcs::write(GFitsHDU* hdu)
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns coordinate system.
 ***************************************************************************/
std::string GWcs::coordsys(void) const
{
    // Set coordinate system
    std::string s_coordsys;
    switch (m_coordsys) {
    case 0:
        s_coordsys = "EQU";
        break;
    case 1:
        s_coordsys = "GAL";
        break;
    default:
        s_coordsys = "UNKNOWN";
        break;
    }

    // Return coordinate system
    return s_coordsys;
}


/***********************************************************************//**
 * @brief Set coordinate system.
 *
 * @param[in] coordsys Coordinate system (EQU/CEL/E/C or GAL/G)
 *
 * @exception GException::wcs_bad_coords 
 *            Invalid coordsys parameter.
 ***************************************************************************/
void GWcs::coordsys(const std::string& coordsys)
{
    // Convert argument to upper case
    std::string ucoordsys = toupper(coordsys);

    // Set coordinate system
    if (ucoordsys == "EQU" || ucoordsys == "CEL" || ucoordsys == "E" ||
        ucoordsys == "C")
        m_coordsys = 0;
    else if (ucoordsys == "GAL" || ucoordsys == "G")
        m_coordsys = 1;
    else
        throw GException::wcs_bad_coords(G_COORDSYS_SET, coordsys);

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                           GWcs private methods                          =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GWcs::init_members(void)
{
    // Initialise members
    m_coordsys =     0;
    m_cdelt[0] =   0.0;
    m_cdelt[1] =   0.0;
    m_crpix[0] =   0.0;
    m_crpix[1] =   0.0;
    m_crval[0] =   0.0;
    m_crval[1] =   0.0;
    m_longpole = 180.0;
    m_latpole  =   0.0;
    m_cd       = GMatrix(2,2);
    m_invcd    = GMatrix(2,2);
    m_pv2      = GVector(21);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] wcs GWcs instance from which members should be copied
 ***************************************************************************/
void GWcs::copy_members(const GWcs& wcs)
{
    // Copy attributes
    m_coordsys = wcs.m_coordsys;
    m_cdelt[0] = wcs.m_cdelt[0];
    m_cdelt[1] = wcs.m_cdelt[1];
    m_crpix[0] = wcs.m_crpix[0];
    m_crpix[1] = wcs.m_crpix[1];
    m_crval[0] = wcs.m_crval[0];
    m_crval[1] = wcs.m_crval[1];
    m_longpole = wcs.m_longpole;
    m_latpole  = wcs.m_latpole;
    m_cd       = wcs.m_cd;
    m_invcd    = wcs.m_invcd;
    m_pv2      = wcs.m_pv2;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GWcs::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               GWcs friends                              =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                       Other functions used by GWcs                      =
 =                                                                         =
 ==========================================================================*/
