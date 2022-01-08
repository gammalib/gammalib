/***************************************************************************
 *            GPulsarEphemerides.cpp - Pulsar ephemerides class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2022 by Juergen Knoedlseder                              *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GPulsarEphemerides.cpp
 * @brief Pulsar ephemerides class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GPulsarEphemerides.hpp"

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
GPulsarEphemerides::GPulsarEphemerides(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] ephemerides Pulsar ephemerides.
 ***************************************************************************/
GPulsarEphemerides::GPulsarEphemerides(const GPulsarEphemerides& ephemerides)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(ephemerides);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GPulsarEphemerides::~GPulsarEphemerides(void)
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
 * @param[in] ephemerides Pulsar ephemerides.
 * @return Pulsar ephemerides.
 ***************************************************************************/
GPulsarEphemerides& GPulsarEphemerides::operator=(const GPulsarEphemerides& ephemerides)
{
    // Execute only if object is not identical
    if (this != &ephemerides) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(ephemerides);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear Pulsar ephemerides
 ***************************************************************************/
void GPulsarEphemerides::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone Pulsar ephemerides
 *
 * @return Pointer to deep copy of Pulsar ephemerides.
 ***************************************************************************/
GPulsarEphemerides* GPulsarEphemerides::clone(void) const
{
    return new GPulsarEphemerides(*this);
}


/***********************************************************************//**
 * @brief Print Pulsar ephemerides
 *
 * @param[in] chatter Chattiness.
 * @return String containing Pulsar ephemerides information.
 *
 * @todo Implement method.
 ***************************************************************************/
std::string GPulsarEphemerides::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GPulsarEphemerides ===");

        // Append information
        // TODO: Add any relevant information

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GPulsarEphemerides::init_members(void)
{
    // Initialise members
    // TODO: Initialise all data members
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] ephemerides Pulsar ephemerides.
 ***************************************************************************/
void GPulsarEphemerides::copy_members(const GPulsarEphemerides& ephemerides)
{
    // Copy members
    // TODO: Copy all data members

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GPulsarEphemerides::free_members(void)
{
    // Return
    return;
}
