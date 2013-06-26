/***************************************************************************
 *        GSkyRegionCircle.cpp - Class implementing a circular sky region       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Michael Mayer                         *
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
 * @file GSkyRegionCircle.hpp
 * @brief Circular sky region implementation
 * @author Michael Mayer
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GSkyRegionCircle.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Prototype __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GSkyRegionCircle::GSkyRegionCircle(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Direction constructor
 *
 * @param[in] centre GSkyDir of Circle center.
 * @param[in] radius Radius of the region [deg].
 ***************************************************************************/
GSkyRegionCircle::GSkyRegionCircle(GSkyDir &centre, const double &radius)
{
    // Initialise members
	init_members();

	// set members
	m_centre = centre;
	m_radius = radius;

	// compute solid angle
	compute_solid();

    // Return
    return;
}

/***********************************************************************//**
 * @brief Direction constructor
 *
 * @param[in] ra Right ascension of Circle center [deg].
 * @param[in] dec declination of Circle center [deg].
 * @param[in] radius Radius of the region [deg].
 ***************************************************************************/
//GSkyRegionCircle::GSkyRegionCircle(const double&  ra, const double&  dec, const double& radius)
//{
//    // Initialise members
//    init_members();
//
//    // Set parameters
//    GSkyDir dir = GSkyDir();
//    dir.radec(ra,dec);
//    m_centre = dir;
//    m_radius = radius;
//
//    // compute solid angle
//    compute_solid();
//
//    // Return
//    return;
//}


GSkyRegionCircle::GSkyRegionCircle(const std::string &line)
{
	 // Initialise members
	 init_members();

	 // read from line
	 read(line);

	 // Return
	 return;

}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] circular sky region .
 ***************************************************************************/
GSkyRegionCircle::GSkyRegionCircle(const GSkyRegionCircle& circle)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(circle);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GSkyRegionCircle::~GSkyRegionCircle(void)
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
 * @param[in] circle SkyRegionCircle object.
 ***************************************************************************/
GSkyRegionCircle& GSkyRegionCircle::operator= (const GSkyRegionCircle& circle)
{
    // Execute only if object is not identical
    if (this != &circle) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(circle);

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
void GSkyRegionCircle::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone object
 ***************************************************************************/
GSkyRegionCircle* GSkyRegionCircle::clone(void) const
{
    // Clone this Region
    return new GSkyRegionCircle(*this);
}


/***********************************************************************//**
 * @brief Set radius value of region
 *
 * @param[in] radius Radius value [deg].
 ***************************************************************************/
void GSkyRegionCircle::radius(const double& radius)
{
    // Set radius value
    m_radius = radius;

    // Recompute the solid angle
    compute_solid();
}


/***********************************************************************//**
 * @brief Set centre of region
 *
 * @param[in] dir Center of region.
 ***************************************************************************/
void GSkyRegionCircle::centre(const GSkyDir& dir)
{
    // Set y value
    m_centre = dir;
}


/***********************************************************************//**
 * @brief Set centre values
 *
 * @param[in] ra Right ascension value.
 * @param[in] dec declination value.
 ***************************************************************************/
void GSkyRegionCircle::centre(const double& ra, const double& dec)
{
    // Set centre values
    m_centre.radec_deg(ra,dec);
}


void GSkyRegionCircle::read(const std::string line) const
{
	std::cout<<"Reading of GSkyRegionCircle::read not implemented yet"<<std::endl;
}
std::string GSkyRegionCircle::write() const
{
	std::cout<<"Writing of GSkyRegionCircle::read not implemented yet"<<std::endl;
	return "";
}
/***********************************************************************//**
 * @brief Print circular region
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing region information.
 ***************************************************************************/
std::string GSkyRegionCircle::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append string
    	result.append("=== GSkyRegionCircle ===");
    	result.append("\n(");
        result.append(gammalib::str(m_centre.ra_deg()));
        result.append(",");
        result.append(gammalib::str(m_centre.dec_deg()));
        result.append(",");
        result.append(gammalib::str(m_radius));
        result.append(")");

    } // endif: chatter was not silent

    // Return result
    return result;
}

/***********************************************************************//**
 * @brief checks if coordinate is within region
 *
 * @param[in] dir Sky direction.
 ***************************************************************************/
bool GSkyRegionCircle::contains(const GSkyDir& dir) const
{
	// Initialise return value
	bool dir_is_in = false;

	// calculate distance from dir to centre
	double distance = dir.dist_deg(m_centre);

	// Check if distance < radius
	if (distance <= radius()) {
		dir_is_in = true;
	}

	// Return bool
	return dir_is_in;
}

/***********************************************************************//**
 * @brief checks if region is fully contained within this region
 *
 * @param[in] reg Sky region.
 ***************************************************************************/
bool GSkyRegionCircle::contains(const GSkyRegion& reg) const
{
	// Initialise return value
	bool fully_inside = false;

	// If other region is circle use a simple way to calculate
	if (reg.type() == "Circle")
	{
		// create circular region from reg
		GSkyRegionCircle* regcirc = (GSkyRegionCircle*)reg.clone();

		// calculate angular distance between the centers
		double ang_dist = m_centre.dist_deg(regcirc->centre());

		// Check if the region is contained in this
		if ((ang_dist + regcirc->radius()) <= m_radius)
		{
			// Set return value to true
			fully_inside = true;
		}
	}
	else
	{
		std::cout<<"ERROR: GSkyRegionCircle::contains(const GSkyRegion& reg) not implemented yet"<<std::endl;
	}
	// Return value
	return fully_inside;
}

/***********************************************************************//**
 * @brief checks if region is overlapping with this region
 *
 * @param[in] reg Sky region.
 ***************************************************************************/
bool GSkyRegionCircle::overlaps(const GSkyRegion& reg) const
{
	// Initialise return value
	bool overlap = false;

	// If other region is circle use a simple way to calculate
	if (reg.type() == "Circle")
	{
		// create circular region from reg
		GSkyRegionCircle* regcirc =  (GSkyRegionCircle*)reg.clone();

		// calculate angular distance between the centers
		double ang_dist = m_centre.dist_deg(regcirc->centre());

		// Check if the distance is smaller than the sum of both radii
		if (ang_dist <= (m_radius + regcirc->radius()))
		{
			// Set overlap to true
			overlap = true;
		}
	}

	else
	{
		std::cout<<"ERROR: GSkyRegionCircle::overlaps(const GSkyRegion& reg) not implemented yet"<<std::endl;
	}

	// Return value
	return overlap;
}



/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GSkyRegionCircle::init_members(void)
{
    // Initialise members
	m_centre = GSkyDir();
	m_radius = 0.0;
	m_type = "Circle";

	//Return
	return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] region Circular sky region.
 ***************************************************************************/
void GSkyRegionCircle::copy_members(const GSkyRegionCircle& region)
{

	// Copy attributes
	m_centre = region.centre();
	m_radius = region.radius();
	m_solid = region.solidangle();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GSkyRegionCircle::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief compute solid angle
 ***************************************************************************/
void GSkyRegionCircle::compute_solid(void)
{
	// convert radius to radians
	double radius_rad = m_radius * gammalib::deg2rad;

	// compute solid angle
	m_solid = gammalib::twopi * (1 - std::cos(radius_rad));

    // Return
    return;
}


