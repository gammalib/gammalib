/***************************************************************************
 *          GLATInstDir.cpp - Fermi/LAT instrument direction class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Jurgen Knodlseder                           *
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
 * @file GLATInstDir.cpp
 * @brief Fermi/LAT instrument direction class implementation
 * @author Juergen Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GLATInstDir.hpp"
#include "GTools.hpp"
#include "GMath.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Prototypes _________________________________________________________ */

/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GLATInstDir::GLATInstDir(void) : GInstDir()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief GSkyDir constructor
 *
 * @param[in] dir Sky direction.
 *
 * Construct LAT instrument direction from sky direction.
 ***************************************************************************/
GLATInstDir::GLATInstDir(const GSkyDir& dir) : GInstDir()
{
    // Initialise class members
    init_members();

    // Assign sky direction
    m_dir = dir;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] dir Instrument direction.
 ***************************************************************************/
GLATInstDir::GLATInstDir(const GLATInstDir& dir) : GInstDir(dir)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(dir);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GLATInstDir::~GLATInstDir(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] dir Instrument direction.
 * @return Instrument direction.
 ***************************************************************************/
GLATInstDir& GLATInstDir::operator=(const GLATInstDir& dir)
{
    // Execute only if object is not identical
    if (this != &dir) {

        // Copy base class members
        this->GInstDir::operator=(dir);

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(dir);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear Fermi/LAT instrument direction
 ***************************************************************************/
void GLATInstDir::clear(void)
{
    // Free members
    free_members();
    this->GInstDir::free_members();

    // Initialise private members
    this->GInstDir::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone Fermi/LAT instrument direction
 *
 * @return Pointer to deep copy of Fermi/LAT instrument direction.
 ***************************************************************************/
GLATInstDir* GLATInstDir::clone(void) const
{
    return new GLATInstDir(*this);
}


/***********************************************************************//**
 * @brief Print instrument direction information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing instrument direction information.
 ***************************************************************************/
std::string GLATInstDir::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append instrument direction
        result.append("RA="+gammalib::str(dir().ra_deg()) +
                      ", DEC="+gammalib::str(dir().dec_deg()));

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
void GLATInstDir::init_members(void)
{
    // Initialise members
    m_dir.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] dir Instrument direction.
 ***************************************************************************/
void GLATInstDir::copy_members(const GLATInstDir& dir)
{
    // Copy attributes
    m_dir = dir.m_dir;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATInstDir::free_members(void)
{
    // Return
    return;
}
