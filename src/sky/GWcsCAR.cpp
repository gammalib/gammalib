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
#include "GException.hpp"
#include "GWcsCAR.hpp"

/* __ Method name definitions ____________________________________________ */

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
 * @param[in] crval Coordinates of reference pixel.
 * @param[in] crpix1 X coordinate of reference pixel.
 * @param[in] crpix2 Y coordinate of reference pixel.
 * @param[in] cdelt1 X increment at reference point in degrees/pixel.
 * @param[in] cdelt2 Y increment at reference point in degrees/pixel.
 * @param[in] coords Coordinate system ('EQU' or 'GAL').
 ***************************************************************************/
GWcsCAR::GWcsCAR(GSkyDir& crval, const double& crpix1, const double& crpix2,
                 const double& cdelt1, const double& cdelt2,
                 const std::string& coords) :
                 GWcs(crval, crpix1, crpix2, cdelt1, cdelt2, coords)
{
    // Initialise class members
    init_members();

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
    os << " Coordinate system .........: " << wcs.coordsys();

    // Return output stream
    return os;
}
