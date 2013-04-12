/***************************************************************************
 *                GCOMPointing.cpp - COMPTEL pointing class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2013 by Juergen Knoedlseder                         *
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
 * @file GCOMPointing.cpp
 * @brief COMPTEL pointing class interface implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GCOMPointing.hpp"

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
GCOMPointing::GCOMPointing(void) : GPointing()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Sky direction constructor
 *
 * @param[in] dir Sky direction.
 *
 * Construct COMPTEL pointing from sky direction.
 ***************************************************************************/
GCOMPointing::GCOMPointing(const GSkyDir& dir) : GPointing()
{
    // Initialise members
    init_members();

    // Assign sky direction
    this->dir(dir);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] pnt COMPTEL pointing.
 ***************************************************************************/
GCOMPointing::GCOMPointing(const GCOMPointing& pnt) : GPointing(pnt)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(pnt);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCOMPointing::~GCOMPointing(void)
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
 * @param[in] pnt COMPTEL pointing.
 ***************************************************************************/
GCOMPointing& GCOMPointing::operator= (const GCOMPointing& pnt)
{
    // Execute only if object is not identical
    if (this != &pnt) {

        // Copy base class members
        this->GPointing::operator=(pnt);

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(pnt);

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
void GCOMPointing::clear(void)
{
    // Free members
    free_members();
    this->GPointing::free_members();

    // Initialise private members
    this->GPointing::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 ***************************************************************************/
GCOMPointing* GCOMPointing::clone(void) const
{
    return new GCOMPointing(*this);
}


/***********************************************************************//**
 * @brief Return pointing direction
 *
 * Returns reference to sky direction of the pointing.
 ***************************************************************************/
const GSkyDir& GCOMPointing::dir(void) const
{
    // Return
    return m_dir;
}


/***********************************************************************//**
 * @brief Set pointing direction
 *
 * @param[in] dir Sky direction.
 *
 * Sets the sky direction for the pointing.
 ***************************************************************************/
void GCOMPointing::dir(const GSkyDir& dir)
{
    // Set sky direction
    m_dir = dir;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print COMPTEL pointing information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing COMPTEL pointing information.
 ***************************************************************************/
std::string GCOMPointing::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMPointing ===");

        // Append information
        result.append("\n"+gammalib::parformat("Pointing direction"));
        result.append(this->dir().print(chatter));

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
void GCOMPointing::init_members(void)
{
    // Initialise members
    m_dir.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] pnt COM pointing.
 ***************************************************************************/
void GCOMPointing::copy_members(const GCOMPointing& pnt)
{
    // Copy members
    m_dir = pnt.m_dir;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMPointing::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
