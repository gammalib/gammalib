/***************************************************************************
 *                GEvent.cpp  -  Event abstract base class                 *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2009 by Jurgen Knodlseder                   *
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
 * @file GEvent.cpp
 * @brief GEvent abstract base class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include <iomanip.h>
#include "GException.hpp"
#include "GEvent.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                       GEvent constructors/destructors                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GEvent::GEvent()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] event Event from which the instance should be built.
 ***************************************************************************/
GEvent::GEvent(const GEvent& event)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(event);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GEvent::~GEvent()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             GEvent operators                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] event Event which should be assigned.
 ***************************************************************************/
GEvent& GEvent::operator= (const GEvent& event)
{
    // Execute only if object is not identical
    if (this != &event) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(event);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                          GEvent public methods                          =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                          GEvent private methods                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GEvent::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] event GEvent members which should be copied.
 ***************************************************************************/
void GEvent::copy_members(const GEvent& event)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GEvent::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                              GEvent friends                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put event into an output stream
 *
 * @param[in] os Output stream into which the event will be dumped
 * @param[in] event Event to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GEvent& event)
{
    // Put event in output stream
    os << (&event)->pipe(os);
        
    // Return output stream
    return os;
}


/*==========================================================================
 =                                                                         =
 =                     Other functions used by GEvent                      =
 =                                                                         =
 ==========================================================================*/
