/***************************************************************************
 *        GSPIInstDir.cpp - INTEGRAL/SPI instrument direction class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GSPIInstDir.cpp
 * @brief INTEGRAL/SPI instrument direction class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstring>  // memcpy
#include "GSPIInstDir.hpp"

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
 *
 * Creates an empty INTEGRAL/SPI instrument direction.
 ***************************************************************************/
GSPIInstDir::GSPIInstDir(void) : GInstDir()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief SPI instrument direction constructor
 *
 * @param[in] dir Pointing direction.
 * @param[in] detid Detector identifier.
 ***************************************************************************/
GSPIInstDir::GSPIInstDir(const GSkyDir& dir, const int& detid) : GInstDir()
{
    // Initialise class members
    init_members();

    // Set members
    m_dir   = dir;
    m_detid = detid;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] dir INTEGRAL/SPI instrument direction.
 ***************************************************************************/
GSPIInstDir::GSPIInstDir(const GSPIInstDir& dir) : GInstDir(dir)
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
GSPIInstDir::~GSPIInstDir(void)
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
 * @param[in] dir INTEGRAL/SPI instrument direction.
 * @return INTEGRAL/SPI instrument direction.
 ***************************************************************************/
GSPIInstDir& GSPIInstDir::operator=(const GSPIInstDir& dir)
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
 * @brief Clear INTEGRAL/SPI instrument direction
 *
 * Clears INTEGRAL/SPI instrument direction by resetting all class members to
 * an initial state. Any information that was present before will be lost.
 ***************************************************************************/
void GSPIInstDir::clear(void)
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
 * @brief Clone INTEGRAL/SPI instrument direction
 *
 * @return Pointer to deep copy of INTEGRAL/SPI instrument direction.
 ***************************************************************************/
GSPIInstDir* GSPIInstDir::clone(void) const
{
    return new GSPIInstDir(*this);
}


/***********************************************************************//**
 * @brief Return instrument direction hash value
 *
 * @return Hash value.
 *
 * Returns a hash value that can be used in the response cache.
 ***************************************************************************/
u_int64_t GSPIInstDir::hash(void) const
{
    // Allocate static array to store the information as floats
    static float buffer[2];

    // Shift detector ID for addition to Right Ascension and Declination in
    // radians
    float shifted_detid = float(m_detid+10);

    // Store the two sky coordinates as floats
    buffer[0] = float(m_dir.ra()  + shifted_detid);
    buffer[1] = float(m_dir.dec() + shifted_detid);

    // Map the floats to an unsigned 64 Bit integer
    u_int64_t hash; std::memcpy(&hash, &buffer, sizeof hash);

    // Return hash value
    return hash;
}


/***********************************************************************//**
 * @brief Print INTEGRAL/SPI instrument direction information
 *
 * @param[in] chatter Chattiness.
 * @return String containing INTEGRAL/SPI instrument direction information.
 ***************************************************************************/
std::string GSPIInstDir::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append instrument direction
        std::string msg = "RA=" + gammalib::str(m_dir.ra_deg()) +
                          ", DEC=" + gammalib::str(m_dir.dec_deg()) +
                          ", DETID=" + gammalib::str(m_detid);
        result.append(msg);

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
void GSPIInstDir::init_members(void)
{
    // Initialise members
    m_dir.clear();
    m_detid = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] dir INTEGRAL/SPI instrument direction.
 ***************************************************************************/
void GSPIInstDir::copy_members(const GSPIInstDir& dir)
{
    // Copy members
    m_dir   = dir.m_dir;
    m_detid = dir.m_detid;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GSPIInstDir::free_members(void)
{
    // Return
    return;
}
