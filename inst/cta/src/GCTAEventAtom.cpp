/***************************************************************************
 *                GCTAEventAtom.cpp - CTA event atom class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2015 by Jurgen Knodlseder                           *
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
 * @file GCTAEventAtom.cpp
 * @brief CTA event atom class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <string>
#include <cmath>
#include "GCTAEventAtom.hpp"
#include "GCTAException.hpp"
#include "GTools.hpp"

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
GCTAEventAtom::GCTAEventAtom(void) : GEventAtom()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] atom Event atom.
 ***************************************************************************/
GCTAEventAtom::GCTAEventAtom(const GCTAEventAtom& atom) : GEventAtom(atom)
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
GCTAEventAtom::~GCTAEventAtom(void)
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
 * @param[in] atom Event atom.
 * @return Event atom.
 ***************************************************************************/
GCTAEventAtom& GCTAEventAtom::operator=(const GCTAEventAtom& atom)
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
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear event atom
 ***************************************************************************/
void GCTAEventAtom::clear(void)
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
 * @return Pointer to deep copy of event atom.
 ***************************************************************************/
GCTAEventAtom* GCTAEventAtom::clone(void) const
{
    return new GCTAEventAtom(*this);
}


/***********************************************************************//**
 * @brief Print event information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing event information.
 ***************************************************************************/
std::string GCTAEventAtom::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append event information
        result.append("Dir="+m_dir.print());
        result.append(" Energy="+m_energy.print());
        result.append(" Time="+m_time.print());

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
void GCTAEventAtom::init_members(void)
{
    // Initialise members
    m_dir.clear();
    m_time.clear();
    m_energy.clear();
    m_index       = 0;
    m_event_id    = 0;
    m_phase       = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] atom Event atom.
 ***************************************************************************/
void GCTAEventAtom::copy_members(const GCTAEventAtom& atom)
{
    // Copy members
    m_dir         = atom.m_dir;
    m_time        = atom.m_time;
    m_energy      = atom.m_energy;
    m_index       = atom.m_index;
    m_event_id    = atom.m_event_id;
    m_phase       = atom.m_phase;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAEventAtom::free_members(void)
{
    // Return
    return;
}
