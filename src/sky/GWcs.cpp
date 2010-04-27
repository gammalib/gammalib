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
 * @brief Constructor
 *
 * @param[in] crval Coordinates of reference pixel.
 * @param[in] crpix1 X coordinate of reference pixel.
 * @param[in] crpix2 Y coordinate of reference pixel.
 * @param[in] cdelt1 X increment at reference point in degrees/pixel.
 * @param[in] cdelt2 Y increment at reference point in degrees/pixel.
 * @param[in] coords Coordinate system ('EQU' or 'GAL').
 ***************************************************************************/
GWcs::GWcs(GSkyDir& crval, const double& crpix1, const double& crpix2,
           const double& cdelt1, const double& cdelt2,
           const std::string& coords)
{
    // Initialise class members
    init_members();

    // Set coordinate system
    coordsys(coords);

    // Set reference value
    switch (m_coordsys) {
    case 0:
        m_crval[0] = crval.ra_deg();
        m_crval[1] = crval.dec_deg();
        break;
    case 1:
        m_crval[0] = crval.l_deg();
        m_crval[1] = crval.b_deg();
        break;
    default:
        //TODO: Throw error
        break;
    }

    // Set members
    m_crpix[0] = crpix1;
    m_crpix[1] = crpix2;
    m_cdelt[0] = cdelt1;
    m_cdelt[1] = cdelt2;

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
 * @exception GException::wcs_bad_coords Invalid coordsys parameter.
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
    m_ctype[0].clear();
    m_ctype[1].clear();
    m_cdelt[0] =   0.0;
    m_cdelt[1] =   0.0;
    m_crpix[0] =   0.0;
    m_crpix[1] =   0.0;
    m_crval[0] =   0.0;
    m_crval[1] =   0.0;
    m_longpole = 180.0;
    m_latpole  =   0.0;
    m_pv2[0]   =   0.0;
    m_pv2[1]   =   0.0;
    m_coordsys =     0;

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
    m_ctype[0] = wcs.m_ctype[0];
    m_ctype[1] = wcs.m_ctype[1];
    m_cdelt[0] = wcs.m_cdelt[0];
    m_cdelt[1] = wcs.m_cdelt[1];
    m_crpix[0] = wcs.m_crpix[0];
    m_crpix[1] = wcs.m_crpix[1];
    m_crval[0] = wcs.m_crval[0];
    m_crval[1] = wcs.m_crval[1];
    m_longpole = wcs.m_longpole;
    m_latpole  = wcs.m_latpole;
    m_pv2[0]   = wcs.m_pv2[0];
    m_pv2[1]   = wcs.m_pv2[1];
    m_coordsys = wcs.m_coordsys;

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
