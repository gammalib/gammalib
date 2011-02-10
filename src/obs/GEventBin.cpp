/***************************************************************************
 *              GEventBin.cpp  -  Event bin abstract base class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GEventBin.cpp
 * @brief Abstract event bin base class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GEventBin.hpp"

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
 * @brief Void constructor
 ***************************************************************************/
GEventBin::GEventBin(void) : GEvent()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] bin Event bin.
 ***************************************************************************/
GEventBin::GEventBin(const GEventBin& bin) : GEvent(bin)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(bin);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GEventBin::~GEventBin(void)
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
 * @param[in] bin Event bin.
 ***************************************************************************/
GEventBin& GEventBin::operator= (const GEventBin& bin)
{
    // Execute only if object is not identical
    if (this != &bin) {

        // Copy base class members
        this->GEvent::operator=(bin);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(bin);

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
 * @brief Return error in number of counts
 *
 * Returns \f$\sqrt(counts)\f$ as the uncertainty in the number of counts
 * in the bin. If counts <= 0 then zero is returned.
 ***************************************************************************/
double GEventBin::error(void) const
{
    // Compute uncertainty
    double error = sqrt(counts());

    // Return error
    return error;
}


/*==========================================================================
 =                                                                         =
 =                              Private methods                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GEventBin::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] bin Event bin.
 ***************************************************************************/
void GEventBin::copy_members(const GEventBin& bin)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GEventBin::free_members(void)
{
    // Return
    return;
}
