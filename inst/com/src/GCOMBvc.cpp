/***************************************************************************
 *         GCOMBvc.cpp - COMPTEL Solar System Barycentre Data class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2022 by Juergen Knodlseder                               *
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
 * @file GCOMBvc.cpp
 * @brief COMPTEL Solar System Barycentre Data class implementation
 * @author Juergen Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GSkyDir.hpp"
#include "GCOMBvc.hpp"

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
GCOMBvc::GCOMBvc(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] bvc COMPTEL Solar System Barycentre Data.
 ***************************************************************************/
GCOMBvc::GCOMBvc(const GCOMBvc& bvc)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(bvc);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCOMBvc::~GCOMBvc(void)
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
 * @param[in] bvc COMPTEL Solar System Barycentre Data.
 * @return COMPTEL Solar System Barycentre Data.
 ***************************************************************************/
GCOMBvc& GCOMBvc::operator=(const GCOMBvc& bvc)
{
    // Execute only if object is not identical
    if (this != &bvc) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(bvc);

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
 * @brief Clear COMPTEL Solar System Barycentre Data
 ***************************************************************************/
void GCOMBvc::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone COMPTEL Solar System Barycentre Data
 *
 * @return Pointer to deep copy of COMPTEL Solar System Barycentre Data.
 ***************************************************************************/
GCOMBvc* GCOMBvc::clone(void) const
{
    return new GCOMBvc(*this);
}


/***********************************************************************//**
 * @brief Return time difference between photon arrival time at CGRO and
 *        the Solar System Barycentre (SSB)
 *
 * @param[in] dir Source position.
 * @return Time difference between photon arrival times (s)
 *
 * Returns the time difference between photon arrival time at CGRO and the
 * Solar System Barycentre (SSB). The arrival time at the SSB is computed
 * by adding the time difference to the photon arrival time as measured by
 * COMPTEL
 *
 * \f[
 *    T_{\rm SSB} = T_{\rm CGRO} + \Delta T
 * \f]
 *
 * The routine implements the algorithm PUL-AL-004 and is inspried from the
 * COMPASS code evpbin02.pulssb.f.
 *
 * It computes
 *
 * \f[
 *    \Delta T = \Delta T_{\rm travel} - \Delta T_{\rm rel} + \Delta T_{\rm BVC}
 * \f]
 *
 * where
 *
 * \f[
 *    \Delta T_{\rm travel} = \left( \vec{SSB} \cdot \vec{n} \right) \times 10^{-6}
 * \f]
 *
 * is the light travel time in seconds between CGRO and the SSB, with
 * \f$\vec{SSB}\f$ being the vector going from the SSB to CGRO, specified
 * using the GCOMBvc::ssb() method, and \f$\vec{n}\f$ is the normalised
 * vector of the source position, provided by the GSkyDir::celvector()
 * method,
 *
 * \f[
 *    \Delta T_{\rm rel} =
 *    -2 R \log \left( 1 + \frac{\Delta T_{\rm travel}}{|\vec{SSB}| * 10^{-6}} \right)
 * \f]
 *
 * is the relativistic delay due to the Sun in seconds, with
 * \f$R=0.49254909 \times 10^{-5}\f$ s, and \f$\Delta T_{\rm BVC}\f$ is the
 * difference in seconds due to the time unit conversion as specified by the
 * GCOMBvc::tdelta() method.
 ***************************************************************************/
double GCOMBvc::tdelta(const GSkyDir& dir) const
{
    // Set constants
    const double radius = 0.49254909e-5;

    // Get celestial vector
    GVector n = dir.celvector();

    // Compute the light travel time from the satellite to SSB along the
    // pulsar direction
    double travt = m_ssb * n * 1.0e-6;

    // Compute the relativistic delay due to the Sun
    double r     = norm(m_ssb) * 1.0e-6;
    double relat = -2.0 * radius * std::log(1.0 + (travt/r));

    // Compute the time difference at SSB
    double tdelta = travt - relat + this->tdelta();

    // Return
    return tdelta;
}


/***********************************************************************//**
 * @brief Print COMPTEL Solar System Barycentre Data
 *
 * @param[in] chatter Chattiness.
 * @return String containing COMPTEL Solar System Barycentre Data
 *         information.
 ***************************************************************************/
std::string GCOMBvc::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMBvc ===");

        // Append information
        result.append("\n"+gammalib::parformat("COMPTEL time"));
        result.append(gammalib::str(m_tjd));
        result.append(":");
        result.append(gammalib::str(m_tics));
        result.append("\n"+gammalib::parformat("MJD"));
        result.append(gammalib::str(m_time.mjd()));
        result.append(" days");
        result.append("\n"+gammalib::parformat("UTC"));
        result.append(m_time.utc());

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
void GCOMBvc::init_members(void)
{
    // Initialise members
    m_time.clear();
    m_ssb.clear();
    m_tjd    = 0;
    m_tics   = 0;
    m_tdelta = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] bvc COMPTEL Solar System Barycentre Data.
 ***************************************************************************/
void GCOMBvc::copy_members(const GCOMBvc& bvc)
{
    // Copy members
    m_time   = bvc.m_time;
    m_tjd    = bvc.m_tjd;
    m_tics   = bvc.m_tics;
    m_ssb    = bvc.m_ssb;
    m_tdelta = bvc.m_tdelta;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMBvc::free_members(void)
{
    // Return
    return;
}
