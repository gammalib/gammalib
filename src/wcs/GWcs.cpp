/***************************************************************************
 *           GWcs.cpp  -  World Coordinate System virtual base class       *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2010 by Jurgen Knodlseder                   *
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
 * @file GWcs.cpp
 * @brief GWcs virtual base class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GWcs.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Local prototypes ___________________________________________________ */

/* __ Constants __________________________________________________________ */

/* __ Static conversion arrays ___________________________________________ */

/*==========================================================================
 =                                                                         =
 =                       GWcs constructors/destructors                     =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GWcs::GWcs()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param wcs GWcs instance which should be used for construction
 ***************************************************************************/
GWcs::GWcs(const GWcs& wcs)
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
GWcs::~GWcs()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                            GWcs operators                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] wcs GWcs instance to be assigned
 ***************************************************************************/
GWcs& GWcs::operator= (const GWcs& wcs)
{
    // Execute only if object is not identical
    if (this != &wcs) {

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
 =                           GWcs public methods                           =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                           GWcs private methods                          =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GWcs::init_members(void)
{
    // Initialise members
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] wcs GWcs instance from which members should be copied
 ***************************************************************************/
void GWcs::copy_members(const GWcs& wcs)
{
    // Copy attributes

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GWcs::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               GWcs friends                              =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                       Other functions used by GWcs                      =
 =                                                                         =
 ==========================================================================*/
