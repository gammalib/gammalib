/***************************************************************************
 *                       GTPLClass.cpp - [WHAT] class                      *
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
 * @file GTPLClass.cpp
 * @brief [WHAT] class implementation
 * @author [AUTHOR]
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTPLClass.hpp"

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
GTPLClass::GTPLClass(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] TPL_OBJECT [WHAT].
 ***************************************************************************/
GTPLClass::GTPLClass(const GTPLClass& TPL_OBJECT)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(TPL_OBJECT);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GTPLClass::~GTPLClass(void)
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
 * @param[in] TPL_OBJECT [WHAT].
 * @return [WHAT].
 ***************************************************************************/
GTPLClass& GTPLClass::operator=(const GTPLClass& TPL_OBJECT)
{
    // Execute only if object is not identical
    if (this != &TPL_OBJECT) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(TPL_OBJECT);

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
 * @brief Clear [WHAT]
 ***************************************************************************/
void GTPLClass::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone [WHAT]
 *
 * @return Pointer to deep copy of [WHAT].
 ***************************************************************************/
GTPLClass* GTPLClass::clone(void) const
{
    return new GTPLClass(*this);
}


/***********************************************************************//**
 * @brief Print [WHAT]
 *
 * @param[in] chatter Chattiness.
 * @return String containing [WHAT] information.
 *
 * @todo Implement method.
 ***************************************************************************/
std::string GTPLClass::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GTPLClass ===");

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
void GTPLClass::init_members(void)
{
    // Initialise members
    // TODO: Initialise all data members
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] TPL_OBJECT [WHAT].
 ***************************************************************************/
void GTPLClass::copy_members(const GTPLClass& TPL_OBJECT)
{
    // Copy members
    // TODO: Copy all data members

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GTPLClass::free_members(void)
{
    // Return
    return;
}
