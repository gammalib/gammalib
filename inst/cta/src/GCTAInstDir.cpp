/***************************************************************************
 *            GCTAInstDir.cpp  -  CTA instrument direction class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GCTAInstDir.cpp
 * @brief GCTAInstDir class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GCTAInstDir.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Prototypes _________________________________________________________ */

/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTAInstDir::GCTAInstDir(void) : GInstDir()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief GSkyDir constructor
 *
 * @param[in] dir Sky direction.
 *
 * Construct CTA instrument direction from sky direction.
 ***************************************************************************/
GCTAInstDir::GCTAInstDir(const GSkyDir& dir) : GInstDir()
{
    // Initialise class members
    init_members();

    // Assign sky direction
    m_dir = dir;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] dir Instrument direction.
 ***************************************************************************/
GCTAInstDir::GCTAInstDir(const GCTAInstDir& dir) : GInstDir(dir)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(dir);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAInstDir::~GCTAInstDir(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] dir Instrument direction.
 ***************************************************************************/
GCTAInstDir& GCTAInstDir::operator= (const GCTAInstDir& dir)
{
    // Execute only if object is not identical
    if (this != &dir) {

        // Copy base class members
        this->GInstDir::operator=(dir);

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(dir);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 ***************************************************************************/
void GCTAInstDir::clear(void)
{
    // Free members
    free_members();
    this->GInstDir::free_members();

    // Initialise private members
    this->GInstDir::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 ***************************************************************************/
GCTAInstDir* GCTAInstDir::clone(void) const
{
    return new GCTAInstDir(*this);
}


/***********************************************************************//**
 * @brief Rotate CTA instrument direction by zenith and azimuth angle
 *
 * @param[in] phi Azimuth angle (deg).
 * @param[in] theta Zenith angle (deg).
 *
 * Rotate CTA instrument direction by a zenith and azimuth angle given in
 * the system of the instrument direction and aligned in celestial
 * coordinates.
 * The azimuth angle is counted counter clockwise from celestial north
 * (this is identical to the astronomical definition of a position angle).
 ***************************************************************************/
void GCTAInstDir::rotate_deg(const double& phi, const double& theta)
{
    // Convert instrument direction into sky direction
    GSkyDir sky = skydir();

    // Rotate sky direction
    sky.rotate_deg(phi, theta);

    // Convert sky direction to instrument direction
    skydir(sky);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute angular distance between instrument directions in radians
 *
 * @param[in] dir Instrument direction.
 ***************************************************************************/
double GCTAInstDir::dist(const GCTAInstDir& dir) const
{
    // Assign sky direction from instrument direction
    GSkyDir sky;
    double  ra  = dir.ra();
    double  dec = dir.dec();
    sky.radec(ra,dec);

    // Compute distance
    double dist = m_dir.dist(sky);

    // Return distance
    return dist;
}


/***********************************************************************//**
 * @brief Compute angular distance between instrument directions in degrees
 *
 * @param[in] dir Instrument direction.
 ***************************************************************************/
double GCTAInstDir::dist_deg(const GCTAInstDir& dir) const
{
    // Return distance in degrees
    return (dist(dir) * rad2deg);
}


/***********************************************************************//**
 * @brief Compute position angle between instrument directions in radians
 *
 * @param[in] dir Instrument direction.
 ***************************************************************************/
double GCTAInstDir::posang(const GCTAInstDir& dir) const
{
    // Assign sky direction from instrument direction
    GSkyDir sky;
    double  ra  = dir.ra();
    double  dec = dir.dec();
    sky.radec(ra,dec);

    // Compute position angle
    double pa = m_dir.posang(sky);

    // Return position angle
    return pa;
}


/***********************************************************************//**
 * @brief Compute position angle between instrument directions in degrees
 *
 * @param[in] dir Instrument direction.
 ***************************************************************************/
double GCTAInstDir::posang_deg(const GCTAInstDir& dir) const
{
    // Return position angle in degrees
    return (posang(dir) * rad2deg);
}


/***********************************************************************//**
 * @brief Print instrument direction information
 ***************************************************************************/
std::string GCTAInstDir::print(void) const
{
    // Initialise result string
    std::string result;

    // Append instrument direction
    result.append("RA="+str(ra_deg())+", DEC="+str(dec_deg()));

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
void GCTAInstDir::init_members(void)
{
    // Initialise members
    m_dir.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] dir Instrument direction.
 ***************************************************************************/
void GCTAInstDir::copy_members(const GCTAInstDir& dir)
{
    // Copy attributes
    m_dir = dir.m_dir;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAInstDir::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
