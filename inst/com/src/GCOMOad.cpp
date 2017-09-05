/***************************************************************************
 *              GCOMOad.cpp - COMPTEL Orbit Aspect Data class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knodlseder                               *
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
 * @file GCOMOad.cpp
 * @brief COMPTEL Orbit Aspect Data class implementation
 * @author Juergen Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GCOMOad.hpp"

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
GCOMOad::GCOMOad(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] oad COMPTEL Orbit Aspect Data.
 ***************************************************************************/
GCOMOad::GCOMOad(const GCOMOad& oad)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(oad);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCOMOad::~GCOMOad(void)
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
 * @param[in] oad COMPTEL Orbit Aspect Data.
 * @return COMPTEL Orbit Aspect Data.
 ***************************************************************************/
GCOMOad& GCOMOad::operator=(const GCOMOad& oad)
{
    // Execute only if object is not identical
    if (this != &oad) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(oad);

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
 * @brief Clear COMPTEL Orbit Aspect Data
 ***************************************************************************/
void GCOMOad::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone COMPTEL Orbit Aspect Data
 *
 * @return Pointer to deep copy of COMPTEL Orbit Aspect Data.
 ***************************************************************************/
GCOMOad* GCOMOad::clone(void) const
{
    return new GCOMOad(*this);
}


/***********************************************************************//**
 * @brief Print COMPTEL Orbit Aspect Data
 *
 * @param[in] chatter Chattiness.
 * @return String containing COMPTEL Orbit Aspect Data information.
 ***************************************************************************/
std::string GCOMOad::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMOad ===");

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
void GCOMOad::init_members(void)
{
    // Initialise members
    m_tstart.clear();
    m_tstop.clear();
    m_tjd  = 0;
    m_tics = 0;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] oad COMPTEL Orbit Aspect Data.
 ***************************************************************************/
void GCOMOad::copy_members(const GCOMOad& oad)
{
    // Copy members
    m_tstart = oad.m_tstart;
    m_tstop  = oad.m_tstop;
    m_tjd    = oad.m_tjd;
    m_tics   = oad.m_tics;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMOad::free_members(void)
{
    // Return
    return;
}
