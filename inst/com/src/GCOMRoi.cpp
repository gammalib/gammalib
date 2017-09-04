/***************************************************************************
 *              GCOMRoi.cpp - COMPTEL region of interest class             *
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
 * @file GCOMRoi.cpp
 * @brief COMPTEL region of interest class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GEvent.hpp"
#include "GCOMRoi.hpp"

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
GCOMRoi::GCOMRoi(void) : GRoi()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] roi COMPTEL region of interest.
 ***************************************************************************/
GCOMRoi::GCOMRoi(const GCOMRoi& roi) : GRoi(roi)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(roi);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Region of interest constructor
 *
 * @param[in] centre Instrument direction centre.
 * @param[in] radius Instrument direction radius.
 * @param[in] phibar_min Phibar minimum (deg).
 * @param[in] phibar_max Phibar maximum (deg).
 ***************************************************************************/
GCOMRoi::GCOMRoi(const GCOMInstDir& centre, const double& radius,
                 const double& phibar_min, const double& phibar_max) : GRoi()

{
    // Initialise class members
    init_members();

    // Set members
    this->centre(centre);
    this->radius(radius);
    this->phibar_min(phibar_min);
    this->phibar_max(phibar_max);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCOMRoi::~GCOMRoi(void)
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
 * @param[in] roi COMPTEL region of interest.
 * @return COMPTEL region of interest.
 ***************************************************************************/
GCOMRoi& GCOMRoi::operator=(const GCOMRoi& roi)
{
    // Execute only if object is not identical
    if (this != &roi) {

        // Copy base class members
        this->GRoi::operator=(roi);

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(roi);

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
 * @brief Clear region of interest
 ***************************************************************************/
void GCOMRoi::clear(void)
{
    // Free members
    free_members();
    this->GRoi::free_members();

    // Initialise private members
    this->GRoi::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone region of interest
 *
 * @return Pointer to deep copy of COMPTEL region of interest.
 ***************************************************************************/
GCOMRoi* GCOMRoi::clone(void) const
{
    return new GCOMRoi(*this);
}


/***********************************************************************//**
 * @brief Check if region of interest contains an event
 *
 * @return True if region of interest contains event, false otherwise.
 ***************************************************************************/
bool GCOMRoi::contains(const GEvent& event) const
{
    // Initialise flag to non-containment
    bool contains = false;

    // Get pointer to COMPTEL instrument direction
    const GCOMInstDir* dir = dynamic_cast<const GCOMInstDir*>(&event.dir());

    // If instrument direction is a COMPTEL instrument direction then check
    // on containment
    if (dir != NULL) {
        if ((m_centre.dir().dist_deg(dir->dir()) <= m_radius) &&
            (m_phibar_min                        <= dir->phibar()) &&
            (m_phibar_max                        >= dir->phibar())) {
            contains = true;
        }
    }

    // Return containment flag
    return contains;
}


/***********************************************************************//**
 * @brief Print region of interest information
 *
 * @param[in] chatter Chattiness.
 * @return String containing region of interest information.
 *
 * @todo Implement method.
 ***************************************************************************/
std::string GCOMRoi::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMRoi ===");

        // Append information
        result.append("\n"+gammalib::parformat("RoI centre"));
        result.append(m_centre.print());
        result.append("\n"+gammalib::parformat("RoI radius"));
        result.append(gammalib::str(m_radius)+" deg");
        result.append("\n"+gammalib::parformat("Phibar range"));
        result.append(gammalib::str(m_phibar_min)+" - ");
        result.append(gammalib::str(m_phibar_max)+" deg");

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
void GCOMRoi::init_members(void)
{
    // Initialise members
    m_centre.clear();
    m_radius     = 0.0;
    m_phibar_min = 0.0;
    m_phibar_max = 0.0;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] roi COMPTEL region of interest.
 ***************************************************************************/
void GCOMRoi::copy_members(const GCOMRoi& roi)
{
    // Copy attributes
    m_centre     = roi.m_centre;
    m_radius     = roi.m_radius;
    m_phibar_min = roi.m_phibar_min;
    m_phibar_max = roi.m_phibar_max;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMRoi::free_members(void)
{
    // Return
    return;
}
