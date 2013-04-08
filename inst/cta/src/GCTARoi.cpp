/***************************************************************************
 *                 GCTARoi.cpp - CTA region of interest class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GCTARoi.cpp
 * @brief GCTARoi class interface implementation.
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GCTARoi.hpp"
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
GCTARoi::GCTARoi(void) : GRoi()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] roi Region of interest.
 ***************************************************************************/
GCTARoi::GCTARoi(const GCTARoi& roi) : GRoi(roi)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(roi);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTARoi::~GCTARoi(void)
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
 * @param[in] roi Region of interest.
 ***************************************************************************/
GCTARoi& GCTARoi::operator= (const GCTARoi& roi)
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
 * @brief Clear instance
 ***************************************************************************/
void GCTARoi::clear(void)
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
 * @brief Clone instance
***************************************************************************/
GCTARoi* GCTARoi::clone(void) const
{
    return new GCTARoi(*this);
}


/***********************************************************************//**
 * @brief Print ROI information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing ROI information.
 ***************************************************************************/
std::string GCTARoi::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTARoi ===");

        // Append information
        result.append("\n"+parformat("ROI centre")+m_centre.print());
        result.append("\n"+parformat("ROI radius")+str(m_radius)+" deg");

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
void GCTARoi::init_members(void)
{
    // Initialise members
    m_centre.clear();
    m_radius = 0.0;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] roi Region of interest from which members should be copied.
 ***************************************************************************/
void GCTARoi::copy_members(const GCTARoi& roi)
{
    // Copy attributes
    m_centre = roi.m_centre;
    m_radius = roi.m_radius;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTARoi::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
