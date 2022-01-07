/***************************************************************************
 *         GCOMBvc.cpp - COMPTEL Solar System Barycentre Data class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2022 by Juergen Knodlseder                               *
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
 * @file GCOMBvc.cpp
 * @brief COMPTEL Solar System Barycentre Data class implementation
 * @author Juergen Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GCOMBvc.hpp"

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
GCOMBvc::GCOMBvc(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] bvc COMPTEL Solar System Barycentre Data.
 ***************************************************************************/
GCOMBvc::GCOMBvc(const GCOMBvc& bvc)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(bvc);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCOMBvc::~GCOMBvc(void)
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
 * @param[in] bvc COMPTEL Solar System Barycentre Data.
 * @return COMPTEL Solar System Barycentre Data.
 ***************************************************************************/
GCOMBvc& GCOMBvc::operator=(const GCOMBvc& bvc)
{
    // Execute only if object is not identical
    if (this != &bvc) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(bvc);

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
 * @brief Clear COMPTEL Solar System Barycentre Data
 ***************************************************************************/
void GCOMBvc::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone COMPTEL Solar System Barycentre Data
 *
 * @return Pointer to deep copy of COMPTEL Solar System Barycentre Data.
 ***************************************************************************/
GCOMBvc* GCOMBvc::clone(void) const
{
    return new GCOMBvc(*this);
}


/***********************************************************************//**
 * @brief Print COMPTEL Solar System Barycentre Data
 *
 * @param[in] chatter Chattiness.
 * @return String containing COMPTEL Solar System Barycentre Data
 *         information.
 ***************************************************************************/
std::string GCOMBvc::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMBvc ===");

        // Append information
        result.append("\n"+gammalib::parformat("COMPTEL time"));
        result.append(gammalib::str(m_tjd));
        result.append(":");
        result.append(gammalib::str(m_tics));
        result.append("\n"+gammalib::parformat("Superpacket MJD range"));
        result.append(gammalib::str(m_tstart.mjd()));
        result.append(" - ");
        result.append(gammalib::str(m_tstop.mjd()));
        result.append(" days");
        result.append("\n"+gammalib::parformat("Superpacket UTC range"));
        result.append(m_tstart.utc());
        result.append(" - ");
        result.append(m_tstop.utc());

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
void GCOMBvc::init_members(void)
{
    // Initialise members
    m_tstart.clear();
    m_tstop.clear();
    m_ssb.clear();
    m_tjd    = 0;
    m_tics   = 0;
    m_tdelta = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] bvc COMPTEL Solar System Barycentre Data.
 ***************************************************************************/
void GCOMBvc::copy_members(const GCOMBvc& bvc)
{
    // Copy members
    m_tstart = bvc.m_tstart;
    m_tstop  = bvc.m_tstop;
    m_tjd    = bvc.m_tjd;
    m_tics   = bvc.m_tics;
    m_ssb    = bvc.m_ssb;
    m_tdelta = bvc.m_tdelta;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMBvc::free_members(void)
{
    // Return
    return;
}
