/***************************************************************************
 *             GCOMSelection.cpp - COMPTEL selection set class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knoedlseder                              *
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
 * @file GCOMSelection.cpp
 * @brief COMPTEL selection set class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GCOMSelection.hpp"
#include "GCOMEventAtom.hpp"

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
GCOMSelection::GCOMSelection(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] select COMPTEL selection set.
 ***************************************************************************/
GCOMSelection::GCOMSelection(const GCOMSelection& select)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(select);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCOMSelection::~GCOMSelection(void)
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
 * @param[in] select COMPTEL selection set.
 * @return COMPTEL selection set.
 ***************************************************************************/
GCOMSelection& GCOMSelection::operator=(const GCOMSelection& select)
{
    // Execute only if object is not identical
    if (this != &select) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(select);

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
 * @brief Clear COMPTEL selection set
 ***************************************************************************/
void GCOMSelection::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone COMPTEL selection set
 *
 * @return Pointer to deep copy of COMPTEL selection set.
 ***************************************************************************/
GCOMSelection* GCOMSelection::clone(void) const
{
    return new GCOMSelection(*this);
}


/***********************************************************************//**
 * @brief Initialise selection statistics
 ***************************************************************************/
void GCOMSelection::init_statistics(void) const
{
    // Initialise statistics members
    m_num_events_checked  = 0;
    m_num_events_used     = 0;
    m_num_events_rejected = 0;
    m_num_e1_min          = 0;
    m_num_e1_max          = 0;
    m_num_e2_min          = 0;
    m_num_e2_max          = 0;
    m_num_tof_min         = 0;
    m_num_tof_max         = 0;
    m_num_psd_min         = 0;
    m_num_psd_max         = 0;
    m_num_zeta_min        = 0;
    m_num_zeta_max        = 0;
    m_num_reflag_min      = 0;
    m_num_reflag_max      = 0;
    m_num_vetoflag_min    = 0;
    m_num_vetoflag_max    = 0;
    m_num_no_scatter      = 0;
    m_num_invalid_modcom  = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Check if event should be used
 *
 * @param[in] event Event.
 * @return True if event should be used, false otherwise.
 ***************************************************************************/
bool GCOMSelection::use_event(const GCOMEventAtom& event) const
{
    // Initialise usage flag
    bool use = true;

    // Increment number of checked events
    m_num_events_checked++;

    // Compute zeta for event
    double zeta = event.eha() - event.phibar();

    // Check for bad minitelescopes
    if (event.modcom() < 1 && event.modcom() > 98) {
        m_num_invalid_modcom++;
        use = false;
    }

    // Check whether the event has a scatter angle determined. This is
    // signaled by a scatter angle -10e20 in radians.
    else if (event.theta() < -1.0e3) {
        m_num_no_scatter++;
        use = false;
    }

    // Apply event selection
    else if (event.e1() < m_e1_min) {
        m_num_e1_min++;
        use = false;
    }
    else if (event.e1() > m_e1_max) {
        m_num_e1_max++;
        use = false;
    }
    else if (event.e2() < m_e2_min) {
        m_num_e2_min++;
        use = false;
    }
    else if (event.e2() > m_e2_max) {
        m_num_e2_max++;
        use = false;
    }
    else if (event.tof() < m_tof_min) {
        m_num_tof_min++;
        use = false;
    }
    else if (event.tof() > m_tof_max) {
        m_num_tof_max++;
        use = false;
    }
    else if (event.psd() < m_psd_min) {
        m_num_psd_min++;
        use = false;
    }
    else if (event.psd() > m_psd_max) {
        m_num_psd_max++;
        use = false;
    }
    else if (zeta < m_zeta_min) {
        m_num_zeta_min++;
        use = false;
    }
    else if (zeta > m_zeta_max) {
        m_num_zeta_max++;
        use = false;
    }
    else if (event.reflag() < m_reflag_min) {
        m_num_reflag_min++;
        use = false;
    }
    else if (event.reflag() > m_reflag_max) {
        m_num_reflag_max++;
        use = false;
    }
    else if (event.veto() < m_vetoflag_min) {
        m_num_vetoflag_min++;
        use = false;
    }
    else if (event.veto() > m_vetoflag_max) {
        m_num_vetoflag_max++;
        use = false;
    }

    // Update acceptance and rejection statistics
    if (use) {
        m_num_events_used++;
    }
    else {
        m_num_events_rejected++;
    }

    // Return usage flag
    return use;
}


/***********************************************************************//**
 * @brief Print COMPTEL selection set
 *
 * @param[in] chatter Chattiness.
 * @return String containing COMPTEL selection set information.
 ***************************************************************************/
std::string GCOMSelection::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMSelection ===");

        // Append number of checked events
        result.append("\n"+gammalib::parformat("Checked events"));
        result.append(gammalib::str(m_num_events_checked));
        result.append("\n"+gammalib::parformat("Accepted events"));
        result.append(gammalib::str(m_num_events_used));
        result.append("\n"+gammalib::parformat("Rejected events"));
        result.append(gammalib::str(m_num_events_rejected));

        // Append E1 selection statistics
        result.append("\n"+gammalib::parformat("E1 selection"));
        result.append(gammalib::str(m_num_e1_min)+" < [");
        result.append(gammalib::str(m_e1_min)+" - ");
        result.append(gammalib::str(m_e1_max)+" MeV] < ");
        result.append(gammalib::str(m_num_e1_max));

        // Append E2 selection statistics
        result.append("\n"+gammalib::parformat("E2 selection"));
        result.append(gammalib::str(m_num_e2_min)+" < [");
        result.append(gammalib::str(m_e2_min)+" - ");
        result.append(gammalib::str(m_e2_max)+" MeV] < ");
        result.append(gammalib::str(m_num_e2_max));

        // Append TOF selection statistics
        result.append("\n"+gammalib::parformat("TOF selection"));
        result.append(gammalib::str(m_num_tof_min)+" < [");
        result.append(gammalib::str(m_tof_min)+" - ");
        result.append(gammalib::str(m_tof_max)+"] < ");
        result.append(gammalib::str(m_num_tof_max));

        // Append PSD selection statistics
        result.append("\n"+gammalib::parformat("PSD selection"));
        result.append(gammalib::str(m_num_psd_min)+" < [");
        result.append(gammalib::str(m_psd_min)+" - ");
        result.append(gammalib::str(m_psd_max)+"] < ");
        result.append(gammalib::str(m_num_psd_max));

        // Append zeta angle selection statistics
        result.append("\n"+gammalib::parformat("Zeta angle selection"));
        result.append(gammalib::str(m_num_zeta_min)+" < [");
        result.append(gammalib::str(m_zeta_min)+" - ");
        result.append(gammalib::str(m_zeta_max)+" deg] < ");
        result.append(gammalib::str(m_num_zeta_max));

        // Append rejection flag selection statistics
        result.append("\n"+gammalib::parformat("Rejection flag selection"));
        result.append(gammalib::str(m_num_reflag_min)+" < [");
        result.append(gammalib::str(m_reflag_min)+" - ");
        result.append(gammalib::str(m_reflag_max)+"] < ");
        result.append(gammalib::str(m_num_reflag_max));

        // Append veto flag selection statistics
        result.append("\n"+gammalib::parformat("Veto flag selection"));
        result.append(gammalib::str(m_num_vetoflag_min)+" < [");
        result.append(gammalib::str(m_vetoflag_min)+" - ");
        result.append(gammalib::str(m_vetoflag_max)+"] < ");
        result.append(gammalib::str(m_num_vetoflag_max));

        // Append other statistics
        result.append("\n"+gammalib::parformat("No scatter angle"));
        result.append(gammalib::str(m_num_no_scatter));
        result.append("\n"+gammalib::parformat("Invalid minitelescope"));
        result.append(gammalib::str(m_num_invalid_modcom));

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
void GCOMSelection::init_members(void)
{
    // Initialise members
    m_e1_min       = 0.070;  //!< Minimum D1 energy deposit (MeV)
    m_e1_max       =  20.0;  //!< Maximum D1 energy deposit (MeV)
    m_e2_min       = 0.650;  //!< Minimum D2 energy deposit (MeV)
    m_e2_max       =  30.0;  //!< Maximum D2 energy deposit (MeV)
    m_tof_min      = 115.0;  //!< Minimum TOF window
    m_tof_max      = 130.0;  //!< Maximum TOF window
    m_psd_min      =   0.0;  //!< Minimum PSD window
    m_psd_max      = 110.0;  //!< Maximum PSD window
    m_zeta_min     =   0.0;  //!< Minimum Earth horizon angle - Phibar window
    m_zeta_max     = 180.0;  //!< Maximum Earth horizon angle - Phibar window
    m_reflag_min   =     1;  //!< Minimum rejection flag
    m_reflag_max   =  1000;  //!< Maximum rejection flag
    m_vetoflag_min =     0;  //!< Minimum veto flag
    m_vetoflag_max =     0;  //!< Maximum veto flag
   
    // Initialise statistics
    init_statistics();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] select COMPTEL selection set.
 ***************************************************************************/
void GCOMSelection::copy_members(const GCOMSelection& select)
{
    // Copy members
    m_e1_min       = select.m_e1_min;
    m_e1_max       = select.m_e1_max;
    m_e2_min       = select.m_e2_min;
    m_e2_max       = select.m_e2_max;
    m_tof_min      = select.m_tof_min;
    m_tof_max      = select.m_tof_max;
    m_psd_min      = select.m_psd_min;
    m_psd_max      = select.m_psd_max;
    m_zeta_min     = select.m_zeta_min;
    m_zeta_max     = select.m_zeta_max;
    m_reflag_min   = select.m_reflag_min;
    m_reflag_max   = select.m_reflag_max;
    m_vetoflag_min = select.m_vetoflag_min;
    m_vetoflag_max = select.m_vetoflag_max;

    // Copy statistics
    m_num_events_checked  = select.m_num_events_checked;
    m_num_events_used     = select.m_num_events_used;
    m_num_events_rejected = select.m_num_events_rejected;
    m_num_e1_min          = select.m_num_e1_min;
    m_num_e1_max          = select.m_num_e1_max;
    m_num_e2_min          = select.m_num_e2_min;
    m_num_e2_max          = select.m_num_e2_max;
    m_num_tof_min         = select.m_num_tof_min;
    m_num_tof_max         = select.m_num_tof_max;
    m_num_psd_min         = select.m_num_psd_min;
    m_num_psd_max         = select.m_num_psd_max;
    m_num_zeta_min        = select.m_num_zeta_min;
    m_num_zeta_max        = select.m_num_zeta_max;
    m_num_reflag_min      = select.m_num_reflag_min;
    m_num_reflag_max      = select.m_num_reflag_max;
    m_num_vetoflag_min    = select.m_num_vetoflag_min;
    m_num_vetoflag_max    = select.m_num_vetoflag_max;
    m_num_no_scatter      = select.m_num_no_scatter;
    m_num_invalid_modcom  = select.m_num_invalid_modcom;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMSelection::free_members(void)
{
    // Return
    return;
}
