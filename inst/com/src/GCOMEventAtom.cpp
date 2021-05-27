/***************************************************************************
 *               GCOMEventAtom.cpp - COMPTEL event atom class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2021 by Juergen Knoedlseder                         *
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
 * @file GCOMEventAtom.cpp
 * @brief COMPTEL event atom class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <string>
#include "GTools.hpp"
#include "GCOMTools.hpp"
#include "GCOMSupport.hpp"
#include "GCOMEventAtom.hpp"

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
 * Creates an empty COMPTEL event atom.
 ***************************************************************************/
GCOMEventAtom::GCOMEventAtom(void) : GEventAtom()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] atom COMPTEL event atom.
 ***************************************************************************/
GCOMEventAtom::GCOMEventAtom(const GCOMEventAtom& atom) : GEventAtom(atom)
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
GCOMEventAtom::~GCOMEventAtom(void)
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
 * @param[in] atom COMPTEL event atom.
 * @return COMPTEL event atom.
 ***************************************************************************/
GCOMEventAtom& GCOMEventAtom::operator=(const GCOMEventAtom& atom)
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
 * Clears COMPTEL event atom by resetting all class members to an
 * initial state. Any information that was present before will be lost.
 ***************************************************************************/
void GCOMEventAtom::clear(void)
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
 * @return Pointer to deep copy of COMPTEL event atom.
 ***************************************************************************/
GCOMEventAtom* GCOMEventAtom::clone(void) const
{
    return new GCOMEventAtom(*this);
}


/***********************************************************************//**
 * @brief Set event time
 *
 * @param[in] tjd Truncated Julian Days (days).
 * @param[in] tics COMPTEL ticks (1/8 ms).
 *
 * Sets the event time from the native COMPTEL time format.
 ***************************************************************************/
void GCOMEventAtom::time(const int& tjd, const int& tics)
{
    // Set event time by converting from the native COMPTEL time format
    m_time = gammalib::com_time(tjd, tics);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print event information
 *
 * @param[in] chatter Chattiness.
 * @return String containing event information.
 ***************************************************************************/
std::string GCOMEventAtom::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append event attributes
        result.append("Chi="+gammalib::str(m_dir.dir().l_deg()));
        result.append(" Psi="+gammalib::str(m_dir.dir().b_deg()));
        result.append(" Phibar="+gammalib::str(m_dir.phibar()));
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
void GCOMEventAtom::init_members(void)
{
    // Initialise members
    m_dir.clear();
    m_energy.clear();
    m_time.clear();
    m_e1     = 0.0;
    m_e2     = 0.0;
    m_phibar = 0.0;
    m_theta  = 0.0;
    m_phi    = 0.0;
    m_eha    = 0.0;
    m_psd    = 0;
    m_tof    = 0;
    m_modcom = 0;
    m_reflag = 0;
    m_veto   = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] atom COMPTEL event atom.
 ***************************************************************************/
void GCOMEventAtom::copy_members(const GCOMEventAtom& atom)
{
    // Copy members
    m_dir    = atom.m_dir;
    m_energy = atom.m_energy;
    m_time   = atom.m_time;
    m_e1     = atom.m_e1;
    m_e2     = atom.m_e2;
    m_phibar = atom.m_phibar;
    m_theta  = atom.m_theta;
    m_phi    = atom.m_phi;
    m_eha    = atom.m_eha;
    m_psd    = atom.m_psd;
    m_tof    = atom.m_tof;
    m_modcom = atom.m_modcom;
    m_reflag = atom.m_reflag;
    m_veto   = atom.m_veto;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMEventAtom::free_members(void)
{
    // Return
    return;
}
