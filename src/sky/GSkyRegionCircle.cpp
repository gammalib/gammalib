/***************************************************************************
 *               GSkyRegionCircle.cpp - Circular sky region class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Michael Mayer                                    *
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
#define G_CONSTRUCTOR "GSkyRegionCircle::GSkyRegionCircle(GSkyDir&, double&)"
#define G_RADIUS                          "GSkyRegionCircle::radius(double&)"
#define G_READ                         "GSkyRegionCircle::read(std::string&)"
#define G_CONTAINS                  "GSkyRegionCircle::contains(GSkyRegion&)"
#define G_OVERLAPS                  "GSkyRegionCircle::overlaps(GSkyRegion&)"

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
 *
 * @exception GException::invalid_argument
 *            Radius value is less than 0.
 ***************************************************************************/
GSkyRegionCircle::GSkyRegionCircle(GSkyDir& centre, const double& radius)
{
    // Initialise members
	init_members();

	// set members
	m_centre = centre;
	m_radius = radius;

	// Check if radius is valid
	if (m_radius < 0.0) {
		throw GException::invalid_argument(G_CONSTRUCTOR,
              "Radius of a region can't be less than 0");
	}

	// compute solid angle
	compute_solid();

    // Return
    return;
}

/***********************************************************************//**
 * @brief string constructor
 *
 * @param[in] line a string containing a line of a ds9 region file
 ***************************************************************************/
GSkyRegionCircle::GSkyRegionCircle(const std::string& line)
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
 * @param[in] circle Circular sky region.
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
 * @param[in] circle Circular sky region.
 * @return Circular sky region.
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
    this->GSkyRegion::free_members();

    // Initialise private members
    this->GSkyRegion::init_members();
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
 *
 * @exception GException::invalid_argument
 *            Radius value is less than 0.
 ***************************************************************************/
void GSkyRegionCircle::radius(const double& radius)
{
    // Set radius value
    m_radius = radius;

	if (m_radius < 0.0) {
		throw GException::invalid_argument(G_RADIUS,
              "Radius of a region can't be less than 0");
	}

    // Recompute the solid angle
    compute_solid();
}


/***********************************************************************//**
 * @brief Read region from string
 *
 * @param[in] line String in ds9 format.
 *
 * @exception GException::invalid_value
 *            Invalid value found in ds9 format string.
 ***************************************************************************/
void GSkyRegionCircle::read(const std::string& line)
{
	// Clear the current instance
	clear();

	// Split the string into 2 parts seperated by #
	std::vector<std::string> substrings = gammalib::split(line,"#");

	std::string region_def = substrings[0];
	std::string comment    = substrings[1];

	// Finding the circle
	if (region_def.find("circle") == std::string::npos) {
		throw GException::invalid_value(G_READ, 
              "Cannnot find the key word \"circle\" in provided string");
	}

	// Get the the coordinate system of the values
	std::string system = gammalib::split(region_def,";")[0];

	// Get the substring of the important values
	unsigned    pos          = region_def.find("circle(");
	unsigned    end          = region_def.find(")");
	std::string circlestring = region_def.substr(pos+7,end);
	circlestring.erase(circlestring.find(")"),1);

	// Get the values of the region x,y,and radius
	std::vector<std::string> values = gammalib::split(circlestring,",");
	if (values.size()!=3) {
		throw GException::invalid_value(G_READ, 
              "Invalid number of arguments, circle takes exactly 3 arguments");
	}
	double x      = gammalib::todouble(values[0]);
	double y      = gammalib::todouble(values[1]);
	double radius = gammalib::todouble(values[2]);

	// Initialise sky dir
	GSkyDir dir = GSkyDir();

	// Set the values in the correct system
	if (system == "fk5") {
		dir.radec_deg(x,y);
	}
	else if (system == "galactic") {
		dir.lb_deg(x,y);
	}
	else {
		std::stringstream s;
		s << "Provided coordinate system \"" << system << "\" is not supported.";
		throw GException::invalid_value(G_READ, s.str());
	}

	// Setting par values
	m_centre = dir;
	m_radius = radius;

	// Throwing exception if radius is less than 0
	if (m_radius < 0.0) {
			throw GException::invalid_value(G_READ,
                  "Radius of a region can't be less than 0");
		}

	// Compute solid angle
	compute_solid();

	// Check if there is a given name for the region and set it
	std::vector<std::string>comments =gammalib::split(comment, " ");
	for (int i = 0; i < comments.size(); i++) {
		if (gammalib::contains(comments[i],"text")){
			std::vector<std::string> attributes = gammalib::split(comments[i],"=");
			if (attributes.size()<2) {
				throw GException::invalid_value(G_READ, 
                      "Invalid number of arguments, type of attribute must"
                      " be key=value");
			}
			m_name = attributes[1];
		}
	}

	// Return
	return;
}


/***********************************************************************//**
 * @brief write region into a string
 *
 * @return String to be written in a ds9 region file
 ***************************************************************************/
std::string GSkyRegionCircle::write(void) const
{
	std::string result;
	result.append("fk5;circle(");
	result.append(gammalib::str(m_centre.ra_deg()));
	result.append(",");
	result.append(gammalib::str(m_centre.dec_deg()));
	result.append(",");
	result.append(gammalib::str(m_radius));
	result.append(") # text=");
	result.append(m_name);
	result.append("\n");
	return result;
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
 *
 * @exception GException::feature_not_implemented
 *            Regions differ in type.
 ***************************************************************************/
bool GSkyRegionCircle::contains(const GSkyRegion& reg) const
{
	// Initialise return value
	bool fully_inside = false;

	// If other region is circle use a simple way to calculate
	if (reg.type() == "Circle") {

		// create circular region from reg
		GSkyRegionCircle* regcirc = (GSkyRegionCircle*)reg.clone();

		// calculate angular distance between the centers
		double ang_dist = m_centre.dist_deg(regcirc->centre());

		// Check if the region is contained in this
		if ((ang_dist + regcirc->radius()) <= m_radius) {
			fully_inside = true;
		}
        
	}
    
    // ... otherwise throw an exception
	else {
		throw GException::feature_not_implemented(G_CONTAINS,
              "Cannot compare two different region types yet");
	}
    
	// Return value
	return fully_inside;
}

/***********************************************************************//**
 * @brief Checks if region is overlapping with this region
 *
 * @param[in] reg Sky region.
 *
 * @exception GException::feature_not_implemented
 *            Regions differ in type.
 ***************************************************************************/
bool GSkyRegionCircle::overlaps(const GSkyRegion& reg) const
{
	// Initialise return value
	bool overlap = false;

	// If other region is circle use a simple way to calculate
	if (reg.type() == "Circle") {

		// Create circular region from reg
		GSkyRegionCircle* regcirc =  (GSkyRegionCircle*)reg.clone();

		// Calculate angular distance between the centers
		double ang_dist = m_centre.dist_deg(regcirc->centre());

		// Check if the distance is smaller than the sum of both radii
		if (ang_dist <= (m_radius + regcirc->radius())) {
			overlap = true;
		}
        
	}

    // ... otherwise throw an exception
	else {
		throw GException::feature_not_implemented(G_OVERLAPS,
              "Cannot compare two different region types yet");
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
	m_type   = "Circle";

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
	m_centre = region.m_centre;
	m_radius = region.m_radius;

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
	m_solid = gammalib::twopi * (1.0 - std::cos(radius_rad));

    // Return
    return;
}
