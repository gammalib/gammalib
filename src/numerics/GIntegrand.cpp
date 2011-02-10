/***************************************************************************
 *       GIntegrand.cpp  -  Integrand function abstract base class         *
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
 * @file GIntegrand.cpp
 * @brief Abstract integrand base class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GIntegrand.hpp"

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
GIntegrand::GIntegrand(void)
{
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] arg Integrand.
 ***************************************************************************/
GIntegrand::GIntegrand(const GIntegrand& arg)
{
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(arg);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GIntegrand::~GIntegrand(void)
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
 * @param[in] arg Integrand.
 ***************************************************************************/
GIntegrand& GIntegrand::operator= (const GIntegrand& arg)
{
    // Execute only if object is not identical
    if (this != &arg) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(arg);

    } // endif: object was not identical

    // Return
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
void GIntegrand::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] arg Integrand.
 ***************************************************************************/
void GIntegrand::copy_members(const GIntegrand& arg)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GIntegrand::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
