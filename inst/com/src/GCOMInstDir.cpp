/***************************************************************************
 *          GCOMInstDir.cpp - COMPTEL instrument direction class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2017 by Juergen Knoedlseder                         *
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
 * @file GCOMInstDir.cpp
 * @brief COMPTEL instrument direction class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstring>  // memcpy
#include "GCOMInstDir.hpp"
#include "GTools.hpp"

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
GCOMInstDir::GCOMInstDir(void) : GInstDir()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] dir Instrument direction.
 ***************************************************************************/
GCOMInstDir::GCOMInstDir(const GCOMInstDir& dir) : GInstDir(dir)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(dir);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Instrument direction constructor
 *
 * @param[in] dir Sky direction.
 * @param[in] phibar Phibar (deg).
 *
 * Constructs a COMPTEL instrument direction from a sky direction and a
 * phibar value.
 ***************************************************************************/
GCOMInstDir::GCOMInstDir(const GSkyDir& dir, const double& phibar) : GInstDir()
{
    // Initialise class members
    init_members();

    // Set class members
    this->dir(dir);
    this->phibar(phibar);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCOMInstDir::~GCOMInstDir(void)
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
GCOMInstDir& GCOMInstDir::operator= (const GCOMInstDir& dir)
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
 * @brief Clear instance
 ***************************************************************************/
void GCOMInstDir::clear(void)
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
 * @brief Clone instance
 *
 * @return Pointer to deep copy of instrument direction.
 ***************************************************************************/
GCOMInstDir* GCOMInstDir::clone(void) const
{
    return new GCOMInstDir(*this);
}


/***********************************************************************//**
 * @brief Return instrument direction hash value
 *
 * @return Hash value.
 *
 * Returns a hash value that can be used in the response cache.
 ***************************************************************************/
u_int64_t GCOMInstDir::hash(void) const
{
    // Allocate static array to store the information as floats
    static float buffer[2];

    // Scale Phibar for addition to Right Ascension and Declination in
    // radians
    float scaled_phibar = 10.0 * float(m_phibar);

    // Store the two sky coordinates as floats
    buffer[0] = float(m_dir.ra()  + scaled_phibar);
    buffer[1] = float(m_dir.dec() + scaled_phibar);

    // Map the floats to an unsigned 64 Bit integer
    u_int64_t hash; std::memcpy(&hash, &buffer, sizeof hash);

    // Return hash value
    return hash;
}


/***********************************************************************//**
 * @brief Print instrument direction information
 *
 * @param[in] chatter Chattiness.
 * @return String containing instrument direction information.
 ***************************************************************************/
std::string GCOMInstDir::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMInstDir ===");

        // Append information
        result.append("\n"+gammalib::parformat("Sky direction (Chi,Psi)"));
        result.append(m_dir.print(chatter));
        result.append("\n"+gammalib::parformat("Scatter angle (Phibar)"));
        result.append(gammalib::str(m_phibar)+" deg");

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
void GCOMInstDir::init_members(void)
{
    // Initialise members
    m_dir.clear();
    m_phibar = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] dir Instrument direction.
 ***************************************************************************/
void GCOMInstDir::copy_members(const GCOMInstDir& dir)
{
    // Copy members
    m_dir    = dir.m_dir;
    m_phibar = dir.m_phibar;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMInstDir::free_members(void)
{
    // Return
    return;
}
