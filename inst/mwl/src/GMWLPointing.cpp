/***************************************************************************
 *          GMWLPointing.cpp  -  Multi-wavelength pointing class           *
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
 * @file GMWLPointing.cpp
 * @brief GMWLPointing class interface implementation.
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GMWLPointing.hpp"
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
GMWLPointing::GMWLPointing(void) : GPointing()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] pnt Pointing.
 ***************************************************************************/
GMWLPointing::GMWLPointing(const GMWLPointing& pnt) : GPointing(pnt)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(pnt);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GMWLPointing::~GMWLPointing(void)
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
 * @param[in] pnt Pointing.
 ***************************************************************************/
GMWLPointing& GMWLPointing::operator= (const GMWLPointing& pnt)
{
    // Execute only if object is not identical
    if (this != &pnt) {

        // Copy base class members
        this->GPointing::operator=(pnt);

        // Free members
        free_members();

        // Initialise private members
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
void GMWLPointing::clear(void)
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
GMWLPointing* GMWLPointing::clone(void) const
{
    return new GMWLPointing(*this);
}


/***********************************************************************//**
 * @brief Print pointing information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing pointing information
 ***************************************************************************/
std::string GMWLPointing::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GMWLPointing ===");

        // Append information
        result.append("\n"+gammalib::parformat("Pointing direction"));
        result.append("not applicable");

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
void GMWLPointing::init_members(void)
{
    // Initialise members
    m_dir.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] pnt Pointing.
 ***************************************************************************/
void GMWLPointing::copy_members(const GMWLPointing& pnt)
{
    // Copy members
    m_dir = pnt.m_dir;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GMWLPointing::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
