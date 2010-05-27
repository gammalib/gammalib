/***************************************************************************
 *            GLATInstDir.cpp  -  LAT instrument direction class           *
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
 * @file GLATInstDir.cpp
 * @brief GLATInstDir class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GLATInstDir.hpp"

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
 * @brief Constructor
 ***************************************************************************/
GLATInstDir::GLATInstDir(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] dir Instrument direction from which class should be
 *                instantiated.
 ***************************************************************************/
GLATInstDir::GLATInstDir(const GLATInstDir& dir)
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
GLATInstDir::~GLATInstDir(void)
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
GLATInstDir& GLATInstDir::operator= (const GLATInstDir& dir)
{
    // Execute only if object is not identical
    if (this != &dir) {

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
 * @brief Compute angular distance between instrument directions in radians
 *
 * @param[in] dir Instrument direction to which distance is to be computed.
 ***************************************************************************/
double GLATInstDir::dist(GLATInstDir& dir) const
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
double GLATInstDir::dist_deg(GLATInstDir& dir) const
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
void GLATInstDir::init_members(void)
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
void GLATInstDir::copy_members(const GLATInstDir& dir)
{
    // Copy attributes
    m_dir = dir.m_dir;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATInstDir::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 ***************************************************************************/
GLATInstDir* GLATInstDir::clone(void) const
{
    return new GLATInstDir(*this);
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
std::ostream& operator<< (std::ostream& os, const GLATInstDir& dir)
{
    // Output sky direction
    os << dir.m_dir;

    // Return output stream
    return os;
}
