/***************************************************************************
 *               GEvents.hpp - Abstract event container class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2013 by Juergen Knoedlseder                         *
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
 * @file GEvents.cpp
 * @brief Abstract event container class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GEvents.hpp"

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
GEvents::GEvents(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] events Event container.
 ***************************************************************************/
GEvents::GEvents(const GEvents& events)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(events);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GEvents::~GEvents(void)
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
 * @param[in] events Event container.
 * @return Event container.
 ***************************************************************************/
GEvents& GEvents::operator=(const GEvents& events)
{
    // Execute only if object is not identical
    if (this != &events) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(events);

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
 * @brief Set energy boundaries
 *
 * @param[in] ebounds Energy boundaries.
 ***************************************************************************/
void GEvents::ebounds(const GEbounds& ebounds)
{
    // Store energy boundaries
    m_ebounds = ebounds;

    // Call (abstract) energy boundary update method
    set_energies();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set Good Time Intervals
 *
 * @param[in] gti Good Time Intervals.
 ***************************************************************************/
void GEvents::gti(const GGti& gti)
{
    // Store Good Time Intervals
    m_gti = gti;

    // Call (abstract) good time interval update method
    set_times();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GEvents::init_members(void)
{
    // Initialise members
    m_ebounds.clear();
    m_gti.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] events Event container.
 ***************************************************************************/
void GEvents::copy_members(const GEvents& events)
{
    // Copy members
    m_ebounds = events.m_ebounds;
    m_gti     = events.m_gti;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GEvents::free_members(void)
{
    // Return
    return;
}
