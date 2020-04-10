/***************************************************************************
 *             GCTAInstDir.cpp - CTA instrument direction class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2020 by Juergen Knoedlseder                         *
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
 * @file GCTAInstDir.cpp
 * @brief CTA instrument direction class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GMath.hpp"
#include "GCTAInstDir.hpp"

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
GCTAInstDir::GCTAInstDir(void) : GInstDir()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Sky direction constructor
 *
 * @param[in] dir Sky direction.
 *
 * Construct CTA instrument direction from sky direction.
 ***************************************************************************/
GCTAInstDir::GCTAInstDir(const GSkyDir& dir) : GInstDir()
{
    // Initialise class members
    init_members();

    // Assign sky direction
    this->dir(dir);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Instrument coordinates constructor
 *
 * @param[in] detx Instrument coordinate X (radians).
 * @param[in] dety Instrument coordinate Y (radians).
 *
 * Construct CTA instrument direction from instrument coordinates.
 ***************************************************************************/
GCTAInstDir::GCTAInstDir(const double& detx, const double& dety) : GInstDir()
{
    // Initialise class members
    init_members();

    // Assign instrument coordinates
    this->detx(detx);
    this->dety(dety);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Instrument direction constructor
 *
 * @param[in] dir Sky direction.
 * @param[in] detx Instrument coordinate X (radians).
 * @param[in] dety Instrument coordinate Y (radians).
 *
 * Construct CTA instrument direction from sky direction and instrument
 * coordinates.
 ***************************************************************************/
GCTAInstDir::GCTAInstDir(const GSkyDir& dir,
                         const double&  detx,
                         const double&  dety) : GInstDir()
{
    // Initialise class members
    init_members();

    // Assign instrument coordinates
    this->dir(dir);
    this->detx(detx);
    this->dety(dety);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] dir CTA instrument direction.
 ***************************************************************************/
GCTAInstDir::GCTAInstDir(const GCTAInstDir& dir) : GInstDir(dir)
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
GCTAInstDir::~GCTAInstDir(void)
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
 * @param[in] dir CTA instrument direction.
 * @return CTA instrument direction.
 ***************************************************************************/
GCTAInstDir& GCTAInstDir::operator=(const GCTAInstDir& dir)
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
 * @brief Clear CTA instrument direction
 ***************************************************************************/
void GCTAInstDir::clear(void)
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
 * @brief CTA instrument direction
 *
 * @return Pointer to deep copy of CTA instrument direction.
 ***************************************************************************/
GCTAInstDir* GCTAInstDir::clone(void) const
{
    return new GCTAInstDir(*this);
}


/***********************************************************************//**
 * @brief Return instrument direction hash value
 *
 * @return Hash value.
 *
 * Returns a hash value that can be used in the response cache.
 ***************************************************************************/
u_int64_t GCTAInstDir::hash(void) const
{
    // Allocate static array to store the information as floats
    static float buffer[2];

    // Store the two sky coordinates as floats
    buffer[0] = float(m_dir.ra());
    buffer[1] = float(m_dir.dec());

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
std::string GCTAInstDir::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append instrument direction
        std::string msg = "RA=" + gammalib::str(m_dir.ra_deg()) +
                          ", DEC=" + gammalib::str(m_dir.dec_deg()) +
                          " [" + gammalib::str(m_detx) +
                          "," + gammalib::str(m_dety) + "]";
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
void GCTAInstDir::init_members(void)
{
    // Initialise members
    m_dir.clear();
    m_detx     = 0.0;
    m_dety     = 0.0;
    m_has_dir  = false;
    m_has_detx = false;
    m_has_dety = false;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] dir CTA instrument direction.
 ***************************************************************************/
void GCTAInstDir::copy_members(const GCTAInstDir& dir)
{
    // Copy attributes
    m_dir      = dir.m_dir;
    m_detx     = dir.m_detx;
    m_dety     = dir.m_dety;
    m_has_dir  = dir.m_has_dir;
    m_has_detx = dir.m_has_detx;
    m_has_dety = dir.m_has_dety;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAInstDir::free_members(void)
{
    // Return
    return;
}
