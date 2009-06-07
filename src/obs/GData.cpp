/***************************************************************************
 *                  GData.cpp  -  Data abstract base class                 *
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
 * @file GData.cpp
 * @brief GData abstract base class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GData.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                       GData constructors/destructors                    =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GData::GData()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param data GData instance which should be used for construction
 ***************************************************************************/
GData::GData(const GData& data)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(data);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GData::~GData()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             GData operators                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] data GData instance to be assigned
 ***************************************************************************/
GData& GData::operator= (const GData& data)
{
    // Execute only if object is not identical
    if (this != &data) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(data);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                           GData public methods                          =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                           GData private methods                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GData::init_members(void)
{
    // Initialise members

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] data GData instance from which members should be copied
 ***************************************************************************/
void GData::copy_members(const GData& data)
{
    // Copy attributes
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GData::free_members(void)
{
    // Free memory

    // Mark memory as free

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               GData friends                             =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                       Other functions used by GData                     =
 =                                                                         =
 ==========================================================================*/
