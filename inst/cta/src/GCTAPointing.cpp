/***************************************************************************
 *                  GCTAPointing.cpp - CTA pointing class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2021 by Juergen Knoedlseder                         *
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
#include "GTools.hpp"
#include "GFilename.hpp"
#include "GMatrix.hpp"
#include "GSkyDir.hpp"
#include "GXmlElement.hpp"
#include "GCTAPointing.hpp"
#include "GCTAInstDir.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ_XML                         "GCTAPointing::read(GXmlElement&)"
#define G_WRITE_XML                       "GCTAPointing::write(GXmlElement&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
//#define G_USE_VECTORS

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
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Construct CTA pointing from XML element.
 ***************************************************************************/
GCTAPointing::GCTAPointing(const GXmlElement& xml)
{
    // Initialise members
    init_members();

    // Read information from XML element
    read(xml);

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

    // Signal that pointing is valid
    m_valid = true;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get instrument direction from sky direction
 *
 * @param[in] skydir Sky direction.
 * @return CTA instrument direction.
 *
 * Returns instrument direction as function of a sky direction
 ***************************************************************************/
GCTAInstDir GCTAPointing::instdir(const GSkyDir& skydir) const
{
    #if defined(G_USE_VECTORS)
	// Compute rotation matrix
	update();

	// Get celestial vector from sky coordinate
	GVector celvector = skydir.celvector();

	// Transform to instrument system
	GVector inst = m_Rback.transpose() * celvector;

	// Initialise instrument coordinates
    double detx(0.0);
    double dety(0.0);

	// Get offset from cel vector, i.e. distance from camera center
	double theta = std::acos(inst[2]);

	// Check if theta and phi are defined
	if (theta > 0.0 ) {
		double phi = std::asin(inst[1] / std::sin(theta));
		detx       = theta * std::cos(phi);
		dety       = theta * std::sin(phi);
	}
    #else
    double theta = m_dir.dist(skydir);
    double phi   = m_dir.posang(skydir); // Celestial system
    double detx(0.0);
    double dety(0.0);
	if (theta > 0.0 ) {
		detx = theta * std::cos(phi);
		dety = theta * std::sin(phi);
	}
    #endif

	// Initialise instrument direction
	GCTAInstDir inst_direction(skydir);
	inst_direction.detx(detx);
	inst_direction.dety(dety);

    // Return
    return inst_direction;
}


/***********************************************************************//**
 * @brief Get sky direction direction from instrument direction
 *
 * @param[in] instdir Instrument direction.
 * @return Sky direction.
 *
 * Returns sky direction as function of an instrument direction
 ***************************************************************************/
GSkyDir GCTAPointing::skydir(const GCTAInstDir& instdir)const
{
	// Compute rotation matrix
	update();

	// Retrieve instrument coordinates
	double inst_x = instdir.detx();
	double inst_y = instdir.dety();

	// Convert to polar coordinates
	double theta     = std::sqrt(inst_x * inst_x + inst_y * inst_y);
	double phi       = std::atan2(inst_y, inst_x);
    double sin_phi   = std::sin(phi);
    double cos_phi   = std::cos(phi);
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);

	// Build vector from polar coordinates
	GVector native(-cos_phi*sin_theta, sin_phi*sin_theta, cos_theta);

	// Rotate from instrument system into sky system
	GVector skyvector = m_Rback * native;
	GSkyDir sky;
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
 * @brief Read pointing from XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::invalid_value
 *            Invalid XML format encountered.
 *
 * Read pointing parameter from an XML element. The format of the pointing
 * parameter is
 *
 *     <parameter name="Pointing" ra="83.0" dec="22.1"/>
 *
 ***************************************************************************/
void GCTAPointing::read(const GXmlElement& xml)
{
    // Clear pointing
    clear();

    // Get pointing parameter
    const GXmlElement* par = gammalib::xml_get_par(G_READ_XML, xml, "Pointing");

    // Extract position attributes
    if (par->has_attribute("ra") && par->has_attribute("dec")) {
        double ra  = gammalib::todouble(par->attribute("ra"));
        double dec = gammalib::todouble(par->attribute("dec"));
        m_dir.radec_deg(ra,dec);
    }
    else if (par->has_attribute("lon") && par->has_attribute("lat")) {
        double lon = gammalib::todouble(par->attribute("lon"));
        double lat = gammalib::todouble(par->attribute("lat"));
        m_dir.lb_deg(lon,lat);
    }
    else {
        std::string msg = "Attributes \"ra\" and \"dec\" or \"lon\" and"
                          " \"lat\"not found in XML parameter \"Pointing\"."
                          " Please verify the XML format.";
        throw GException::invalid_value(G_READ_XML, msg);
    }

    // Signal that pointing is valid
    m_valid = true;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write pointing information into XML element
 *
 * @param[in] xml XML element.
 *
 * Writes pointing parameter into an XML element. The format of the pointing
 * parameter is
 *
 *     <parameter name="Pointing" ra="83.0" dec="22.1"/>
 *
 ***************************************************************************/
void GCTAPointing::write(GXmlElement& xml) const
{
    // Get parameter
    GXmlElement* par = gammalib::xml_need_par(G_WRITE_XML, xml, "Pointing");

    // Write attributes
    par->attribute("ra",  gammalib::str(m_dir.ra_deg()));
    par->attribute("dec", gammalib::str(m_dir.dec_deg()));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print CTA pointing information
 *
 * @param[in] chatter Chattiness.
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
    m_valid   = false;
    m_zenith  = 0.0;
    m_azimuth = 0.0;

    // Initialise cache
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
    m_valid     = pnt.m_valid;
    m_zenith    = pnt.m_zenith;
    m_azimuth   = pnt.m_azimuth;

    // Copy cache
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
