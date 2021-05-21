/***************************************************************************
 *                GLATEventAtom.cpp - LAT event atom class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2021 by Juergen Knoedlseder                         *
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
 * @file GLATEventAtom.cpp
 * @brief GLATEventAtom class implementation.
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <string>
#include "GTools.hpp"
#include "GLATEventAtom.hpp"

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
 ***************************************************************************/
GLATEventAtom::GLATEventAtom(void) : GEventAtom()
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
GLATEventAtom::GLATEventAtom(const GLATEventAtom& atom) : GEventAtom(atom)
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
GLATEventAtom::~GLATEventAtom(void)
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
GLATEventAtom& GLATEventAtom::operator=(const GLATEventAtom& atom)
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
 ***************************************************************************/
void GLATEventAtom::clear(void)
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
 * @return Pointer to deep copy of event atom
 ***************************************************************************/
GLATEventAtom* GLATEventAtom::clone(void) const
{
    return new GLATEventAtom(*this);
}


/***********************************************************************//**
 * @brief Print event information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing event information.
 ***************************************************************************/
std::string GLATEventAtom::print(const GChatter& chatter) const
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
void GLATEventAtom::init_members(void)
{
    // Initialise members
    m_dir.clear();
    m_time.clear();
    m_energy.clear();
    m_theta               = 0.0;
    m_phi                 = 0.0;
    m_zenith_angle        = 0.0;
    m_earth_azimuth_angle = 0.0;
    m_event_id            = 0;
    m_run_id              = 0;
    m_recon_version       = 0;
    m_calib_version[0]    = 0;
    m_calib_version[1]    = 0;
    m_calib_version[2]    = 0;
    m_event_class         = 0;
    m_conversion_type     = 0;
    m_livetime            = 0.0;
    m_difrsp              = NULL;
    m_num_difrsp          = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] atom Event atom.
 ***************************************************************************/
void GLATEventAtom::copy_members(const GLATEventAtom& atom)
{
    // Copy members
    m_dir                 = atom.m_dir;
    m_time                = atom.m_time;
    m_energy              = atom.m_energy;
    m_theta               = atom.m_theta;
    m_phi                 = atom.m_phi;
    m_zenith_angle        = atom.m_zenith_angle;
    m_earth_azimuth_angle = atom.m_earth_azimuth_angle;
    m_event_id            = atom.m_event_id;
    m_run_id              = atom.m_run_id;
    m_recon_version       = atom.m_recon_version;
    m_calib_version[0]    = atom.m_calib_version[0];
    m_calib_version[1]    = atom.m_calib_version[1];
    m_calib_version[2]    = atom.m_calib_version[2];
    m_event_class         = atom.m_event_class;
    m_conversion_type     = atom.m_conversion_type;
    m_livetime            = atom.m_livetime;

    // Copy other attributes
    m_num_difrsp = atom.m_num_difrsp;

    // If there are diffuse response components then copy them
    if (m_num_difrsp > 0 && atom.m_difrsp != NULL) {

        // Allocate memory for diffuse response components
        m_difrsp = new double[m_num_difrsp];

        // Copy diffuse response components
        for (int i = 0; i < m_num_difrsp; ++i) {
            m_difrsp[i] = atom.m_difrsp[i];
        }

    } // endif: there were diffuse response components to copy

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATEventAtom::free_members(void)
{
    // Free memory
    if (m_difrsp != NULL) delete [] m_difrsp;

    // Signal free pointers
    m_difrsp = NULL;

    // Return
    return;
}
