/***************************************************************************
 *                  GCTAPointing.cpp - CTA pointing class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
#include "GFits.hpp"
#include "GHorizDir.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_DIR_HORIZ                         "GCTAPointing::dir_horiz(GTime&)"

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
 * @brief Pointing table constructor
 *
 * @param[in] filename Pointing table file.
 * @param[in] extname Pointing extension name (defaults to "POINTING")
 *
 * Construct CTA pointing from a pointing table file.
 ***************************************************************************/
GCTAPointing::GCTAPointing(const std::string& filename,
                           const std::string& extname)
{
    // Initialise members
    init_members();

    // Load pointing table
    load(filename, extname);

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
    double phi   = m_dir.posang(skydir);
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


/************************************************************************
 * @brief Return horizontal direction as function of time
 *
 * @param[in] time Time.
 * @return Horizontal direction.
 *
 * @exception GException::out_of_range
 *            Specified time is not in valid range.
 *
 * Returns the horizontal direction for the specified @p time.
 ************************************************************************/
GHorizDir GCTAPointing::dir_horiz(const GTime& time) const
{
    // Initialize horizontal direction
    GHorizDir dir;

    // If no table has been loaded then ... TBD
    if (m_has_table == false) {
        // no pointing table, so either throw exception or return the
        // average direction...
    }

    // ... otherwise interpolate the position
    else {

        // First check if time is inside table bounds

        if (time < m_table_tmin || time > m_table_tmax) {
            throw GException::out_of_range(G_DIR_HORIZ, time.secs(), 
                                           m_table_tmin.secs(),
                                           m_table_tmax.secs());
        }

        // Get interpolated alt and az for the given time using the
        // pointing table
        double alt = m_table_nodes.interpolate(time.secs(), m_table_alt);
        double az  = m_table_nodes.interpolate(time.secs(), m_table_az);
  
        // Set direction
        dir.altaz(alt,az);
    }

    // Return direction
    return dir;
}


/************************************************************************
 * @brief Load pointing table from a pointing file
 *
 * @param[in] filename Pointing table filename.
 * @param[in] extname Pointing extension name (defaults to "POINTING")
 *
 * Loads a pointing table from a FITS file. See the read() method for
 * more information about the structure of the FITS file.
 ************************************************************************/
void GCTAPointing::load(const std::string& filename, const std::string& extname)
{
    // Allocate FITS file
    GFits file;

    // Open FITS file
    file.open(filename);

    // Get pointing table
    const GFitsTable& table = *file.table(extname);

    // Read pointing from table
    read(table);

    // Close FITS file
    file.close();

    // Return
    return;
}


/************************************************************************
 * @brief Load pointing table from a pointing file
 *
 * @param[in] filename Pointing table filename.
 * @param[in] extname Pointing extension name (defaults to "POINTING")
 *
 * Opens a FITS table with columns: START, STOP, ALT_PNT, AZ_PNT,
 * describing the pointing direction as a function of time in
 * horizontal coordinates, and creates an interpolation function for
 * getting the pointing alt/az for an arbitrary time.
 * 
 * The advantage of this method is that ctools does not need to
 * implement ra/dec to alt/az coordinate conversions, since the values
 * are pre-calculated.
 ************************************************************************/
void GCTAPointing::read(const GFitsTable& table)
{
    // Read time reference
    m_reference.read(table);

    // Get number of elements in table
    int nrows = table.nrows();

    // Clear lookup table
    m_table_nodes.clear();
    m_table_az.clear();
    m_table_alt.clear();
    m_table_nodes.reserve(nrows);
    m_table_az.reserve(nrows);
    m_table_alt.reserve(nrows);

    // Get relevant columns
    const GFitsTableCol *start = table["START"];
    const GFitsTableCol *stop  = table["STOP"];
    const GFitsTableCol *alt   = table["ALT_PNT"];
    const GFitsTableCol *az    = table["AZ_PNT"];

    // later on, may also interpolate RA/Dec, in the case of a drift
    //  scan or other tracking mode:
    //     GFitsTableCol *ra     = (*table)["RA_PNT"];
    //     GFitsTableCol *dec    = (*table)["DEC_PNT"];

    // Loop over the elements and build a lookup table
    for (int i = 0; i < nrows; ++i) {

        // Set time
        double midtime = 0.5*(stop->real(i) + start->real(i));
        GTime  time(midtime, m_reference);
        
        // Append to lookup table
        m_table_nodes.append(time.secs());
        m_table_az.push_back(az->real(i)  * gammalib::deg2rad);
        m_table_alt.push_back(alt->real(i) * gammalib::deg2rad);

        // Set minimum and maximum time
        if (i == 0) {
            m_table_tmin = time;
        }
        else if (i == nrows-1) {
            m_table_tmax = time;
        }

    } // endfor: looped over table

    // Signal that table is available
    m_has_table = true;

    // Return
    return;
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
    m_has_table = false;
    m_table_nodes.clear();
    m_table_az.clear();
    m_table_alt.clear();
    m_table_tmin.clear();
    m_table_tmax.clear();
    m_reference.clear();

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
    m_dir         = pnt.m_dir;
    m_zenith      = pnt.m_zenith;
    m_azimuth     = pnt.m_azimuth;
    m_has_table   = pnt.m_has_table;
    m_table_nodes = pnt.m_table_nodes;
    m_table_az    = pnt.m_table_az;
    m_table_alt   = pnt.m_table_alt;
    m_table_tmin  = pnt.m_table_tmin;
    m_table_tmax  = pnt.m_table_tmax;
    m_reference   = pnt.m_reference;

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
