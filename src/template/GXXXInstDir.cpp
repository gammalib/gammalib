/***************************************************************************
 *        GXXXInstDir.cpp - [INSTRUMENT] instrument direction class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) [YEAR] by [AUTHOR]                                       *
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
 * @file GXXXInstDir.cpp
 * @brief [INSTRUMENT] instrument direction class implementation
 * @author [AUTHOR]
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GXXXInstDir.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Prototypes _________________________________________________________ */

/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Creates an empty [INSTRUMENT] instrument direction.
 ***************************************************************************/
GXXXInstDir::GXXXInstDir(void) : GInstDir()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] dir [INSTRUMENT] instrument direction.
 ***************************************************************************/
GXXXInstDir::GXXXInstDir(const GXXXInstDir& dir) : GInstDir(dir)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(dir);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GXXXInstDir::~GXXXInstDir(void)
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
 * @param[in] dir [INSTRUMENT] instrument direction.
 * @return [INSTRUMENT] instrument direction.
 ***************************************************************************/
GXXXInstDir& GXXXInstDir::operator=(const GXXXInstDir& dir)
{
    // Execute only if object is not identical
    if (this != &dir) {

        // Copy base class members
        this->GInstDir::operator=(dir);

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(dir);

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
 * @brief Clear [INSTRUMENT] instrument direction
 *
 * Clears [INSTRUMENT] instrument direction by resetting all class members to
 * an initial state. Any information that was present before will be lost.
 ***************************************************************************/
void GXXXInstDir::clear(void)
{
    // Free members
    free_members();
    this->GInstDir::free_members();

    // Initialise private members
    this->GInstDir::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone [INSTRUMENT] instrument direction
 *
 * @return Pointer to deep copy of [INSTRUMENT] instrument direction.
 ***************************************************************************/
GXXXInstDir* GXXXInstDir::clone(void) const
{
    return new GXXXInstDir(*this);
}


/***********************************************************************//**
 * @brief Print [INSTRUMENT] instrument direction information
 *
 * @param[in] chatter Chattiness.
 * @return String containing [INSTRUMENT] instrument direction information.
 ***************************************************************************/
std::string GXXXInstDir::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append instrument direction
        // TODO: Add code to append any information that might be relevant.
        // Keep the information short, and don't print a header, since you
        // may have many events

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
void GXXXInstDir::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] dir [INSTRUMENT] instrument direction.
 ***************************************************************************/
void GXXXInstDir::copy_members(const GXXXInstDir& dir)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GXXXInstDir::free_members(void)
{
    // Return
    return;
}
