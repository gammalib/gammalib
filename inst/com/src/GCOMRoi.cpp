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
 *
 * @todo Implement method.
 ***************************************************************************/
bool GCOMRoi::contains(const GEvent& event) const
{
    // Initialise flag to non-containment
    bool contains = false;

    // TODO: Implement containment test

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
        // TODO: Add any relevant information

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
    // TODO: Initialise all data members
    // Example:
    double m_radius = 0.0;
    
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
    // TODO: Copy all data members
    // Example:
    m_radius = roi.m_radius;

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
