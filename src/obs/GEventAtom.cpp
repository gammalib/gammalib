/***************************************************************************
 *             GEventAtom.cpp  -  Event atom abstract base class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GEvent.cpp
 * @brief GEvent abstract base class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include <iomanip>
#include "GException.hpp"
#include "GEventAtom.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GEventAtom::GEventAtom(void) : GEvent()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] atom Event atom from which the instance should be built.
 ***************************************************************************/
GEventAtom::GEventAtom(const GEventAtom& atom) : GEvent(atom)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(atom);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GEventAtom::~GEventAtom(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] atom Event atom to be assigned.
 ***************************************************************************/
GEventAtom& GEventAtom::operator= (const GEventAtom& atom)
{
    // Execute only if object is not identical
    if (this != &atom) {

        // Copy base class members
        this->GEvent::operator=(atom);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(atom);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GEventAtom::init_members(void)
{
    // Initialise attributes
    m_time.mjd(0.0);
    m_energy.MeV(0.0);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] atom GEventAtom members which should be copied.
 ***************************************************************************/
void GEventAtom::copy_members(const GEventAtom& atom)
{
    // Copy attributes
    m_time   = atom.m_time;
    m_energy = atom.m_energy;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GEventAtom::free_members(void)
{
    // Return
    return;
}
