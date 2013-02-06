/***************************************************************************
 *                 GXXXPointing.cpp  -  XXX pointing class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
 * @file GXXXPointing.cpp
 * @brief XXX pointing class interface implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GXXXPointing.hpp"

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
GXXXPointing::GXXXPointing(void) : GPointing()
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
 * Construct XXX pointing from sky direction.
 ***************************************************************************/
GXXXPointing::GXXXPointing(const GSkyDir& dir) : GPointing()
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
 * @param[in] pnt XXX pointing.
 ***************************************************************************/
GXXXPointing::GXXXPointing(const GXXXPointing& pnt) : GPointing(pnt)
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
GXXXPointing::~GXXXPointing(void)
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
 * @param[in] pnt XXX pointing.
 * @return XXX pointing.
 ***************************************************************************/
GXXXPointing& GXXXPointing::operator= (const GXXXPointing& pnt)
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
void GXXXPointing::clear(void)
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
 *
 * @return Deep copy of XXX pointing.
 ***************************************************************************/
GXXXPointing* GXXXPointing::clone(void) const
{
    return new GXXXPointing(*this);
}


/***********************************************************************//**
 * @brief Return pointing direction
 *
 * @return Reference to pointing sky direction.
 *
 * Returns reference to sky direction of the pointing.
 ***************************************************************************/
const GSkyDir& GXXXPointing::dir(void) const
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
void GXXXPointing::dir(const GSkyDir& dir)
{
    // Set sky direction
    m_dir = dir;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print XXX pointing information
 * @return String containing XXX pointing information.
 ***************************************************************************/
std::string GXXXPointing::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GXXXPointing ===");

    // Append pointing information
    result.append("\n"+parformat("Pointing direction")+this->dir().print());

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
void GXXXPointing::init_members(void)
{
    // Initialise members
    m_dir.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] pnt XXX pointing.
 ***************************************************************************/
void GXXXPointing::copy_members(const GXXXPointing& pnt)
{
    // Copy members
    m_dir = pnt.m_dir;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GXXXPointing::free_members(void)
{
    // Return
    return;
}
