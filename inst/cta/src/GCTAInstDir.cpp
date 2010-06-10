/***************************************************************************
 *            GCTAInstDir.cpp  -  CTA instrument direction class           *
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
 * @param[in] dir Sky direction from which object is to be constructed.
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
 * @param[in] dir Instrument direction from which class should be
 *                instantiated.
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
 * @param[in] dir Instrument direction to be assigned.
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
 * @brief Clear sky direction
 ***************************************************************************/
void GCTAInstDir::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute angular distance between instrument directions in radians
 *
 * @param[in] dir Instrument direction to which distance is to be computed.
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
 * @param[in] dir Instrument direction to which distance is to be computed.
 ***************************************************************************/
double GCTAInstDir::dist_deg(const GCTAInstDir& dir) const
{
    // Return distance in degrees
    return (dist(dir) * rad2deg);
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
 * @param[in] dir Instrument direction from which members should be copied
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


/***********************************************************************//**
 * @brief Clone instance
 ***************************************************************************/
GCTAInstDir* GCTAInstDir::clone(void) const
{
    return new GCTAInstDir(*this);
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream
 * @param[in] column Instrument direction to put in output stream
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GCTAInstDir& dir)
{
    // Output sky direction
    os << dir.m_dir;

    // Return output stream
    return os;
}
