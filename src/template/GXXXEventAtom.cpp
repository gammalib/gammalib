/***************************************************************************
 *             GXXXEventAtom.cpp - [INSTRUMENT] event atom class           *
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
 * @file GXXXEventAtom.cpp
 * @brief [INSTRUMENT] event atom class implementation
 * @author [AUTHOR]
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <string>
#include "GXXXEventAtom.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Creates an empty [INSTRUMENT] event atom.
 ***************************************************************************/
GXXXEventAtom::GXXXEventAtom(void) : GEventAtom()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] atom [INSTRUMENT] event atom.
 ***************************************************************************/
GXXXEventAtom::GXXXEventAtom(const GXXXEventAtom& atom) : GEventAtom(atom)
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
GXXXEventAtom::~GXXXEventAtom(void)
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
 * @param[in] atom [INSTRUMENT] event atom.
 * @return [INSTRUMENT] event atom.
 ***************************************************************************/
GXXXEventAtom& GXXXEventAtom::operator=(const GXXXEventAtom& atom)
{
    // Execute only if object is not identical
    if (this != &atom) {

        // Copy base class members
        this->GEventAtom::operator=(atom);

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
 =                            Public methods                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear event atom
 *
 * Clears [INSTRUMENT] event atom by resetting all class members to an
 * initial state. Any information that was present before will be lost.
 ***************************************************************************/
void GXXXEventAtom::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GEventAtom::free_members();
    this->GEvent::free_members();

    // Initialise members
    this->GEvent::init_members();
    this->GEventAtom::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone event atom
 *
 * @return Pointer to deep copy of [INSTRUMENT] event atom.
 ***************************************************************************/
GXXXEventAtom* GXXXEventAtom::clone(void) const
{
    return new GXXXEventAtom(*this);
}


/***********************************************************************//**
 * @brief Print event information
 *
 * @param[in] chatter Chattiness.
 * @return String containing event information.
 ***************************************************************************/
std::string GXXXEventAtom::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append event attributes
        result.append("Dir="+m_dir.print(chatter));
        result.append(" Energy="+m_energy.print(chatter));
        result.append(" Time="+m_time.print(chatter));

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GXXXEventAtom::init_members(void)
{
    // Initialise members
    m_dir.clear();
    m_time.clear();
    m_energy.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] atom [INSTRUMENT] event atom.
 ***************************************************************************/
void GXXXEventAtom::copy_members(const GXXXEventAtom& atom)
{
    // Copy members
    m_dir    = atom.m_dir;
    m_time   = atom.m_time;
    m_energy = atom.m_energy;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GXXXEventAtom::free_members(void)
{
    // Return
    return;
}
