/***************************************************************************
 *               GEvents.cpp  -  Events abstract base class                *
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
 * @file GEvents.cpp
 * @brief GEvents abstract base class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GException.hpp"
#include "GEvents.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                      GEvents constructors/destructors                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GEvents::GEvents()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] events Events from which the instance should be built.
 ***************************************************************************/
GEvents::GEvents(const GEvents& events)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(events);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GEvents::~GEvents()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                            GEvents operators                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] events Events which should be assigned.
 ***************************************************************************/
GEvents& GEvents::operator= (const GEvents& events)
{
    // Execute only if object is not identical
    if (this != &events) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(events);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                         GEvents public methods                          =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                         GEvents private methods                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GEvents::init_members(void)
{
    // Initialise members

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] events GEvents members which should be copied.
 ***************************************************************************/
void GEvents::copy_members(const GEvents& events)
{
    // Copy attributes

    // Clone members

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GEvents::free_members(void)
{
    // Free memory

    // Signal free pointers

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             GEvents friends                             =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                     Other functions used by GEvents                     =
 =                                                                         =
 ==========================================================================*/
