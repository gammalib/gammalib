/***************************************************************************
 *                  GCTAPointing.cpp - CTA pointing class                  *
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
 * @file GCTAPointing.cpp
 * @brief CTA pointing class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GCTAPointing.hpp"
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
GCTAPointing::GCTAPointing(void)
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
 * Construct CTA pointing from sky direction.
 ***************************************************************************/
GCTAPointing::GCTAPointing(const GSkyDir& dir)
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
 * @param[in] pnt CTA pointing.
 ***************************************************************************/
GCTAPointing::GCTAPointing(const GCTAPointing& pnt)
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
GCTAPointing::~GCTAPointing(void)
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
 * @param[in] pnt CTA pointing.
 * @return CTA pointing.
 ***************************************************************************/
GCTAPointing& GCTAPointing::operator=(const GCTAPointing& pnt)
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
 * @brief Clear CTA pointing
 ***************************************************************************/
void GCTAPointing::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone CTA pointing
 *
 * @return Poiter to deep copy of CTA pointing.
 ***************************************************************************/
GCTAPointing* GCTAPointing::clone(void) const
{
    return new GCTAPointing(*this);
}


/***********************************************************************//**
 * @brief Set pointing direction
 *
 * @param[in] dir Sky direction of pointing.
 *
 * Set the pointing direction to the specified @p sky direction.
 ***************************************************************************/
void GCTAPointing::dir(const GSkyDir& dir)
{
    // Set sky direction
    m_dir = dir;

    // Invalidate cache
    m_has_cache = false;

    // Return
    return;
}

/***********************************************************************//**
 * @brief Get instrument direction
 *
 * @param[in] dir Sky direction.
 *
 * Returns instrument direction as function of a sky direction
 ***************************************************************************/
const GCTAInstDir& GCTAPointing::instdir(const GSkyDir& skydir) const
{

	// Get celestial vector from sky coordinate
	GVector celvector = skydir.celvector();

	// Transform to instrument system
	GVector inst = m_Rback.transpose() * celvector;

	// Initialise instrument coordinates
	double x_inst = 0.0;
	double y_inst = 0.0;

	// Get offset from cel vector, i.e. distance from camera center
	double theta = std::acos(inst[2]);

	// Check if theta and phi are defined
	if ( theta > 0.0 ) {

		double phi = std::asin( inst[1] / std::sin(theta) );
		x_inst = theta * std::cos(phi);
		y_inst = theta * std::sin(phi);
	}

	// Initialise instrument direction
	GCTAInstDir inst_direction = GCTAInstDir(skydir);
	inst_direction.detx(x_inst);
	inst_direction.dety(y_inst);

    // Return
    return inst_direction;
}


/***********************************************************************//**
 * @brief Get instrument direction
 *
 * @param[in] dir Sky direction.
 *
 * Returns instrument direction as function of a sky direction
 ***************************************************************************/
const GSkyDir& GCTAPointing::skydir(const GCTAInstDir& instdir) const
{

	double inst_x = instdir.detx();
	double inst_y = instdir.dety();

	double theta = std::sqrt( inst_x * inst_x + inst_y * inst_y );
	double phi = std::atan(inst_y / inst_x);

	// Get celestial vector from sky coordinate
	// TODO: this has to be checked!
	// Do we retrieve the correct vector here?
	GVector  native(std::cos(phi)*std::cos(theta),std::cos(phi) * std::sin(theta), std::sin(phi));

	// Transform to instrument system
	GVector skyvector = m_Rback * native;

	GSkyDir sky = GSkyDir();
	sky.celvector(skyvector);

	// Return
	return sky;

}




/***********************************************************************//**
 * @brief Return rotation matrix
 *
 * @return Rotation matrix.
 ***************************************************************************/
const GMatrix& GCTAPointing::rot(void) const
{
    // Update cache
    update();

    // Return rotation matrix
    return m_Rback;
}


/***********************************************************************//**
 * @brief Print CTA pointing information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing pointing information.
 ***************************************************************************/
std::string GCTAPointing::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAPointing ===");

        // Append information
        result.append("\n"+gammalib::parformat("Pointing direction"));
        result.append(this->dir().print());

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
void GCTAPointing::init_members(void)
{
    // Initialise members
    m_dir.clear();
    m_zenith    = 0.0;
    m_azimuth   = 0.0;
    m_has_cache = false;
    m_Rback.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] pnt CTA pointing.
 ***************************************************************************/
void GCTAPointing::copy_members(const GCTAPointing& pnt)
{
    // Copy members
    m_dir       = pnt.m_dir;
    m_zenith    = pnt.m_zenith;
    m_azimuth   = pnt.m_azimuth;
    m_has_cache = pnt.m_has_cache;
    m_Rback     = pnt.m_Rback;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAPointing::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update coordinate transformation cache
 ***************************************************************************/
void GCTAPointing::update(void) const
{
    // Update coordinate transformation cache only if required
    if (!m_has_cache) {

        // Set up Euler matrices
        GMatrix Ry;
        GMatrix Rz;
        Ry.eulery(m_dir.dec_deg() - 90.0);
        Rz.eulerz(-m_dir.ra_deg());

        // Compute rotation matrix
        m_Rback = (Ry * Rz).transpose();

        // Signal that we have a valid transformation cache
        m_has_cache = true;

    } // endif: Update of cache was required

    // Return
    return;
}
