/***************************************************************************
 *                  GCOMDri.cpp - COMPTEL Data Space class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knoedlseder                              *
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
 * @file GCOMDri.cpp
 * @brief COMPTEL Data Space class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GCOMDri.hpp"

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
GCOMDri::GCOMDri(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] dri COMPTEL Data Space.
 ***************************************************************************/
GCOMDri::GCOMDri(const GCOMDri& dri)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(dri);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCOMDri::~GCOMDri(void)
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
 * @param[in] dri COMPTEL Data Space.
 * @return COMPTEL Data Space.
 ***************************************************************************/
GCOMDri& GCOMDri::operator=(const GCOMDri& dri)
{
    // Execute only if object is not identical
    if (this != &dri) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(dri);

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
 * @brief Clear COMPTEL Data Space
 ***************************************************************************/
void GCOMDri::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone COMPTEL Data Space
 *
 * @return Pointer to deep copy of COMPTEL Data Space.
 ***************************************************************************/
GCOMDri* GCOMDri::clone(void) const
{
    return new GCOMDri(*this);
}


/***********************************************************************//**
 * @brief Print COMPTEL Data Space
 *
 * @param[in] chatter Chattiness.
 * @return String containing COMPTEL Data Space information.
 *
 * @todo Implement method.
 ***************************************************************************/
std::string GCOMDri::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMDri ===");

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
void GCOMDri::init_members(void)
{
    // Initialise members
    // TODO: Initialise all data members
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] dri COMPTEL Data Space.
 ***************************************************************************/
void GCOMDri::copy_members(const GCOMDri& dri)
{
    // Copy members
    // TODO: Copy all data members

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMDri::free_members(void)
{
    // Return
    return;
}
