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
 * @param[in] centre Centre sky direction.
 * @param[in] radius Region radius [deg].
 ***************************************************************************/
GSkyRegionCircle::GSkyRegionCircle(GSkyDir& centre, const double& radius)
{
    // Initialise members
	init_members();

	// Set members
	this->centre(centre);
    this->radius(radius);

	// Compute solid angle
	compute_solid_angle();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Direction constructor
 *
 * @param[in] ra Right Ascension of region centre [deg].
 * @param[in] dec Declination of region centre [deg].
 * @param[in] radius Region radius [deg].
 ***************************************************************************/
GSkyRegionCircle::GSkyRegionCircle(const double& ra, const double& dec,
                                   const double& radius)
{
    // Initialise members
	init_members();

	// Set members
	this->centre(ra, dec);
    this->radius(radius);

	// Compute solid angle
	compute_solid_angle();

    // Return
    return;
}


/***********************************************************************//**
 * @brief String constructor
 *
 * @param[in] line DS9 region file line.
 *
 * Constructs region from a DS9 region file line.
 ***************************************************************************/
GSkyRegionCircle::GSkyRegionCircle(const std::string& line)
{
	 // Initialise members
	 init_members();

	 // Read information from DS9 region file line
	 read(line);

	 // Return
	 return;

}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] circle Circular sky region.
 ***************************************************************************/
GSkyRegionCircle::GSkyRegionCircle(const GSkyRegionCircle& region)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(region);

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
GSkyRegionCircle& GSkyRegionCircle::operator=(const GSkyRegionCircle& circle)
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
 * @brief Clone circular sky region
 *
 * @return Deep copy to circular sky region.
 ***************************************************************************/
GSkyRegionCircle* GSkyRegionCircle::clone(void) const
{
    // Clone circular sky region
    return new GSkyRegionCircle(*this);
}


/***********************************************************************//**
 * @brief Set radius of circular region
 *
 * @param[in] radius Radius [deg].
 *
 * @exception GException::invalid_argument
 *            Radius value is less than 0.
 *
 * Sets the radius of the circular sky region. Only non-negative radii are
 * allowed.
 ***************************************************************************/
void GSkyRegionCircle::radius(const double& radius)
{
	if (radius < 0.0) {
        std::string msg =
            "A radius of "+gammalib::str(radius)+" degrees has been"
            " specified for a circular sky region.\n"
            "The radius of a circular sky region can't be less than 0"
            " degrees.";
		throw GException::invalid_argument(G_RADIUS, msg);
	}

    // Set radius value
    m_radius = radius;

    // Compute the solid angle
    compute_solid_angle();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read region from DS9 string
 *
 * @param[in] line String in DS9 format.
 *
 * @exception GException::invalid_value
 *            Invalid value found in DS9 format string.
 ***************************************************************************/
void GSkyRegionCircle::read(const std::string& line)
{
	// Clear the current instance
	clear();

	// Split the string into 2 parts seperated by #
	std::vector<std::string> substrings = gammalib::split(line,"#");
	std::string region_def = (substrings.size() > 0) ? substrings[0] : "";
	std::string comment    = (substrings.size() > 1) ? substrings[1] : "";

	// Finding the circle
	if (region_def.find("circle") == std::string::npos) {
        std::string msg =
            "Unable to find the key word \"circle\" in provided string"
            " \""+line+"\".\n"
            "The \"circle\" key word is mandatory.";
		throw GException::invalid_value(G_READ, msg);
	}

	// Get the the coordinate system of the values
	std::string system = gammalib::split(region_def, ";")[0];

	// Get the substring of the important values
	unsigned    pos          = region_def.find("circle(");
	unsigned    end          = region_def.find(")");
	std::string circlestring = region_def.substr(pos+7, end);
	circlestring.erase(circlestring.find(")"), 1);

	// Get the values of the region x,y,and radius
	std::vector<std::string> values = gammalib::split(circlestring,",");
	if (values.size() != 3) {
        std::string msg =
            "Invalid number of "+gammalib::str(values.size())+" arguments"
            " after the \"circle\" key word in provided string \""+line+"\".\n"
            "Exactly 3 arguments are expected.";
		throw GException::invalid_value(G_READ, msg);
	}
	double x      = gammalib::todouble(values[0]);
	double y      = gammalib::todouble(values[1]);
	double radius = gammalib::todouble(values[2]);

	// Initialise centre direction
	GSkyDir centre = GSkyDir();

	// Set the values in the correct system
	if (system == "fk5") {
		centre.radec_deg(x,y);
	}
	else if (system == "galactic") {
		centre.lb_deg(x,y);
	}
	else {
        std::string msg =
            "Unsupported coordinate system \""+system+"\" in provided string"
            " \""+line+"\".\n"
            "Only the following coordinate systems are supported: \"fk5\", "
            "\"galactic\".";
		throw GException::invalid_value(G_READ, msg);
	}

	// Set members
	this->centre(centre);
    this->radius(radius);

	// Compute solid angle
	compute_solid_angle();

	// Check if there is a given name for the region and set it
	std::vector<std::string>comments = gammalib::split(comment, " ");
	for (int i = 0; i < comments.size(); i++) {
		if (gammalib::contains(comments[i], "text")) {
			std::vector<std::string> attributes = gammalib::split(comments[i], "=");
			if (attributes.size() < 2) {
                std::string msg =
                      "Invalid character sequence encountered in provided"
                      " string \""+line+"\".\n"
                      "An attribute of the type \"text=Name\" is expected.";
				throw GException::invalid_value(G_READ, msg);
			}
			m_name = attributes[1];
		}
	}

	// Return
	return;
}


/***********************************************************************//**
 * @brief Write region into a string
 *
 * @return String to be written in a DS9 region file
 *
 * Writes a DS9 region into a string. The region name is only written if it
 * is defined.
 ***************************************************************************/
std::string GSkyRegionCircle::write(void) const
{
    // Allocate string
	std::string result;

    // Set string
	result.append("fk5;circle(");
	result.append(gammalib::str(m_centre.ra_deg()));
	result.append(",");
	result.append(gammalib::str(m_centre.dec_deg()));
	result.append(",");
	result.append(gammalib::str(m_radius));
	result.append(")");

    // Optionally add region name
    if (m_name.length() > 0) {
        result.append(" # text=");
        result.append(m_name);
    }

    // Return string
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
 * @brief Checks if sky direction lies within region
 *
 * @param[in] dir Sky direction.
 *
 * A sky direction lies within a region when its distance to the region
 * centre is not larger than the region radius.
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
 * @brief Checks if region is fully contained within this region
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

		// Create circular region from reg
		const GSkyRegionCircle* regcirc =
              dynamic_cast<const GSkyRegionCircle*>(&reg);

		// Calculate angular distance between the centers
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
		const GSkyRegionCircle* regcirc =
              dynamic_cast<const GSkyRegionCircle*>(&reg);

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
 * @brief Compute solid angle
 ***************************************************************************/
void GSkyRegionCircle::compute_solid_angle(void)
{
	// Convert radius from degrees to radians
	double radius_rad = m_radius * gammalib::deg2rad;

	// Compute solid angle
	m_solid = gammalib::twopi * (1.0 - std::cos(radius_rad));

    // Return
    return;
}
