/***************************************************************************
 *             GHealpix.cpp  -  Healpix sky representation class           *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
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
 * @file GHealpix.cpp
 * @brief GHealpix class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GHealpix.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                     GHealpix constructors/destructors                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GHealpix::GHealpix()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param pixels GHealpix instance which should be used for construction
 ***************************************************************************/
GHealpix::GHealpix(const GHealpix& pixels)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(pixels);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GHealpix::~GHealpix()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                          GHealpix operators                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] pixels GHealpix instance to be assigned
 ***************************************************************************/
GHealpix& GHealpix::operator= (const GHealpix& pixels)
{
    // Execute only if object is not identical
    if (this != &pixels) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(pixels);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                         GHealpix public methods                         =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                         GHealpix private methods                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GHealpix::init_members(void)
{
    // Initialise members

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] pixels GHealpix instance from which members should be copied
 ***************************************************************************/
void GHealpix::copy_members(const GHealpix& pixels)
{
    // Copy attributes
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GHealpix::free_members(void)
{
    // Free memory

    // Mark memory as free

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone Healpix representation
 ***************************************************************************/
GHealpix* GHealpix::clone(void) const
{
    return new GHealpix(*this);
}


/*==========================================================================
 =                                                                         =
 =                             GHealpix friends                            =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                     Other functions used by GHealpix                    =
 =                                                                         =
 ==========================================================================*/
