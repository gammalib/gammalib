/***************************************************************************
 *                        GPhoton.hpp - Photon class                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
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
 * @file GPhoton.cpp
 * @brief Photon class implementation.
 * @author Juergen Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GPhoton.hpp"
#include "GTools.hpp"

/* __ Constants __________________________________________________________ */

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
GPhoton::GPhoton(void)
{
    // Initialise private members
    init_members();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Photon constructor
 *
 * @param[in] dir Sky direction.
 * @param[in] energy Energy.
 * @param[in] time Time.
 ***************************************************************************/
GPhoton::GPhoton(const GSkyDir& dir, const GEnergy& energy, const GTime& time)
{ 
    // Initialise private members
    init_members();

    // Set members
    m_dir    = dir;
    m_energy = energy;
    m_time   = time;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] ph Photon.
 ***************************************************************************/
GPhoton::GPhoton(const GPhoton& ph)
{ 
    // Initialise private members
    init_members();

    // Copy members
    copy_members(ph);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GPhoton::~GPhoton(void)
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
 * @param[in] ph Photon.
 ***************************************************************************/
GPhoton& GPhoton::operator= (const GPhoton& ph)
{ 
    // Execute only if object is not identical
    if (this != &ph) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(ph);

    } // endif: object was not identical
  
    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 ***************************************************************************/
void GPhoton::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();
    
    // Return
    return; 
}


/***********************************************************************//**
 * @brief Clone object
 ***************************************************************************/
GPhoton* GPhoton::clone(void) const
{
    // Clone this image
    return new GPhoton(*this);
}


/***********************************************************************//**
 * @brief Print photon
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing photon information.
 ***************************************************************************/
std::string GPhoton::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Build photon string
        result.append("GPhoton(");
        result.append("RA="+str(m_dir.ra_deg()));
        result.append(", Dec="+str(m_dir.dec_deg()));
        result.append(", E="+m_energy.print());
        result.append(", MET="+m_time.print());
        if (m_mc_id >= 0) {
            result.append(", MCid="+str(m_mc_id));
        }
        result.append(")");

    } // endif: chatter was not silent

    // Return
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
void GPhoton::init_members(void)
{
    // Initialise members
    m_dir.clear();
    m_energy.clear();
    m_time.clear();
    m_mc_id = -1;
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] ph Photon.
 ***************************************************************************/
void GPhoton::copy_members(const GPhoton& ph)
{
    // Copy time
    m_dir    = ph.m_dir;
    m_energy = ph.m_energy;
    m_time   = ph.m_time;
    m_mc_id  = ph.m_mc_id;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GPhoton::free_members(void)
{
    // Return
    return;
}
