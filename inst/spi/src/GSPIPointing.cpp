/***************************************************************************
 *                 GSPIPointing.cpp  -  SPI pointing class                 *
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
 * @file GSPIPointing.cpp
 * @brief SPI pointing class interface implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GSPIPointing.hpp"

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
GSPIPointing::GSPIPointing(void) : GPointing()
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
 * Construct SPI pointing from sky direction.
 ***************************************************************************/
GSPIPointing::GSPIPointing(const GSkyDir& dir) : GPointing()
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
 * @param[in] pnt SPI pointing.
 ***************************************************************************/
GSPIPointing::GSPIPointing(const GSPIPointing& pnt) : GPointing(pnt)
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
GSPIPointing::~GSPIPointing(void)
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
 * @param[in] pnt SPI pointing.
 * @return SPI pointing.
 ***************************************************************************/
GSPIPointing& GSPIPointing::operator= (const GSPIPointing& pnt)
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
void GSPIPointing::clear(void)
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
 * @return Deep copy of SPI pointing.
 ***************************************************************************/
GSPIPointing* GSPIPointing::clone(void) const
{
    return new GSPIPointing(*this);
}


/***********************************************************************//**
 * @brief Return pointing direction
 *
 * @return Reference to pointing sky direction.
 *
 * Returns reference to sky direction of the pointing.
 ***************************************************************************/
const GSkyDir& GSPIPointing::dir(void) const
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
void GSPIPointing::dir(const GSkyDir& dir)
{
    // Set sky direction
    m_dir = dir;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print SPI pointing information
 * @return String containing SPI pointing information.
 ***************************************************************************/
std::string GSPIPointing::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GSPIPointing ===");
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
void GSPIPointing::init_members(void)
{
    // Initialise members
    m_dir.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] pnt SPI pointing.
 ***************************************************************************/
void GSPIPointing::copy_members(const GSPIPointing& pnt)
{
    // Copy members
    m_dir = pnt.m_dir;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GSPIPointing::free_members(void)
{
    // Return
    return;
}
