/***************************************************************************
 *                 GLATPointing.cpp  -  LAT pointing class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @file GLATPointing.cpp
 * @brief GLATPointing class interface implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GLATPointing.hpp"
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
 * @brief Constructor
 ***************************************************************************/
GLATPointing::GLATPointing(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] pnt Pointing.
 ***************************************************************************/
GLATPointing::GLATPointing(const GLATPointing& pnt)
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
GLATPointing::~GLATPointing(void)
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
GLATPointing& GLATPointing::operator= (const GLATPointing& pnt)
{
    // Execute only if object is not identical
    if (this != &pnt) {

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
void GLATPointing::clear(void)
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
 * @brief Clone class
***************************************************************************/
GLATPointing* GLATPointing::clone(void) const
{
    return new GLATPointing(*this);
}


/***********************************************************************//**
 * @brief Print LAT pointing information
 ***************************************************************************/
std::string GLATPointing::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GLATPointing ===");
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
void GLATPointing::init_members(void)
{
    // Initialise members
    m_dir.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] pnt Pointing from which members should be copied.
 ***************************************************************************/
void GLATPointing::copy_members(const GLATPointing& pnt)
{
    // Copy members
    m_dir = pnt.m_dir;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATPointing::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
