/***************************************************************************
 *                   GEphemerides.cpp - Ephemerides class                  *
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
 * @file GEphemerides.cpp
 * @brief Ephemerides class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GEphemerides.hpp"

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
GEphemerides::GEphemerides(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] ephemerides Ephemerides.
 ***************************************************************************/
GEphemerides::GEphemerides(const GEphemerides& ephemerides)
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
GEphemerides::~GEphemerides(void)
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
 * @param[in] ephemerides Ephemerides.
 * @return Ephemerides.
 ***************************************************************************/
GEphemerides& GEphemerides::operator=(const GEphemerides& ephemerides)
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
 * @brief Clear Ephemerides
 ***************************************************************************/
void GEphemerides::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone Ephemerides
 *
 * @return Pointer to deep copy of Ephemerides.
 ***************************************************************************/
GEphemerides* GEphemerides::clone(void) const
{
    return new GEphemerides(*this);
}


/***********************************************************************//**
 * @brief Print Ephemerides
 *
 * @param[in] chatter Chattiness.
 * @return String containing Ephemerides information.
 ***************************************************************************/
std::string GEphemerides::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GEphemerides ===");

        // Append information
        result.append("\n"+gammalib::parformat("Ephemerides name"));
        result.append(m_name);
        result.append("\n"+gammalib::parformat("File name"));
        result.append(m_filename.url());

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
void GEphemerides::init_members(void)
{
    // Initialise members
    m_name.clear();
    m_filename.clear();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] ephemerides Ephemerides.
 ***************************************************************************/
void GEphemerides::copy_members(const GEphemerides& ephemerides)
{
    // Copy members
    m_name     = ephemerides.m_name;
    m_filename = ephemerides.m_filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GEphemerides::free_members(void)
{
    // Return
    return;
}
