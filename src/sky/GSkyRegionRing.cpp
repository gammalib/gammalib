/***************************************************************************
 *               GSkyRegionRing.cpp - Ring sky region class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Maria Krause, Anneli Schulz                      *
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
 * @file GSkyRegionRing.hpp
 * @brief Ring sky region implementation
 * @author Maria Krause, Anneli Schulz
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GSkyRegionRing.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_RADIUS                          "GSkyRegionRing::radius(double&)"
#define G_READ                         "GSkyRegionRing::read(std::string&)"
#define G_CONTAINS                  "GSkyRegionRing::contains(GSkyRegion&)"
#define G_OVERLAPS                  "GSkyRegionRing::overlaps(GSkyRegion&)"

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
GSkyRegionRing::GSkyRegionRing(void)
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
GSkyRegionRing::GSkyRegionRing(GSkyDir& centre, const double& radius1, const double& radius2)
{
    // Initialise members
	init_members();

	// Set members
	this->centre(centre);
    this->radius1(radius1);
    this->radius2(radius2);

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
 * @param[in] radius1 inner Region radius [deg].
 * @param[in] radius2 outer Region radius [deg].
 ***************************************************************************/
GSkyRegionRing::GSkyRegionRing(const double& ra, const double& dec,
                                   const double& radius1, const double& radius2)
{
    // Initialise members
	init_members();

	// Set members
	this->centre(ra, dec);
    this->radius1(radius1);
    this->radius2(radius2);

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
GSkyRegionRing::GSkyRegionRing(const std::string& line)
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
 * @param[in] region ring sky region.
 ***************************************************************************/
GSkyRegionRing::GSkyRegionRing(const GSkyRegionRing& region)
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
GSkyRegionRing::~GSkyRegionRing(void)
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
 * @param[in] Ring Circular sky region.
 * @return Circular sky region.
 ***************************************************************************/
GSkyRegionRing& GSkyRegionRing::operator=(const GSkyRegionRing& Ring)
{
    // Execute only if object is not identical
    if (this != &Ring) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(Ring);

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
void GSkyRegionRing::clear(void)
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
GSkyRegionRing* GSkyRegionRing::clone(void) const
{
    // Clone circular sky region
    return new GSkyRegionRing(*this);
}


/***********************************************************************//**
 * @brief Set radius of circular region
 *
 * @param[in] radius1 inner region Radius [deg].
 * @param[in] radius2 outer region Radius [deg]. 
 *
 * @exception GException::invalid_argument
 *            Radius value is less than 0.
 *
 * Sets the radius of the circular sky region. Only non-negative radii are
 * allowed.
 ***************************************************************************/
void GSkyRegionRing::radius1(const double& radius1)
{
	if (radius1 < 0.0) {
        std::string msg =
            "A radius of "+gammalib::str(radius1)+" degrees has been"
            " specified for a ring sky region.\n"
            "The radius of a ring sky region can't be less than 0"
            " degrees.";
		throw GException::invalid_argument(G_RADIUS, msg);
	}

    // Set radius value
    m_radius1 = radius1;

    // Compute the solid angle
    compute_solid_angle();

    // Return
    return;
}

void GSkyRegionRing::radius2(const double& radius2)
{
	if (radius2 < 0.0 && radius2 < radius1) {
        std::string msg =
            "A radius of "+gammalib::str(radius2)+" degrees has been"
            " specified for a ring sky region.\n"
            "The radius of a ring sky region can't be less than 0"
            " degrees.";
		throw GException::invalid_argument(G_RADIUS, msg);
	}

    // Set radius value
    m_radius2 = radius2;

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
void GSkyRegionRing::read(const std::string& line)
{
	// Clear the current instance
	clear();

	// Split the string into 2 parts seperated by #
	std::vector<std::string> substrings = gammalib::split(line,"#");
	std::string region_def = (substrings.size() > 0) ? substrings[0] : "";
	std::string comment    = (substrings.size() > 1) ? substrings[1] : "";

	// Finding the Ring
	if (region_def.find("Ring") == std::string::npos) {
        std::string msg =
            "Unable to find the key word \"Ring\" in provided string"
            " \""+line+"\".\n"
            "The \"Ring\" key word is mandatory.";
		throw GException::invalid_value(G_READ, msg);
	}

	// Get the the coordinate system of the values
	std::string system = gammalib::split(region_def, ";")[0];

	// Get the substring of the important values
	unsigned    pos          = region_def.find("Ring(");
	unsigned    end          = region_def.find(")");
	std::string Ringstring = region_def.substr(pos+7, end);
	Ringstring.erase(Ringstring.find(")"), 1);

	// Get the values of the region x,y,and radius
	std::vector<std::string> values = gammalib::split(Ringstring,",");
	if (values.size() != 4) {
        std::string msg =
            "Invalid number of "+gammalib::str(values.size())+" arguments"
            " after the \"Ring\" key word in provided string \""+line+"\".\n"
            "Exactly 4 arguments are expected.";
		throw GException::invalid_value(G_READ, msg);
	}
	double x       = gammalib::todouble(values[0]);
	double y       = gammalib::todouble(values[1]);
	double radius1 = gammalib::todouble(values[2]);
	double radius2 = gammalib::todouble(values[3]);

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
    this->radius1(radius1);
    this->radius2(radius2);

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
std::string GSkyRegionRing::write(void) const
{
    // Allocate string
	std::string result;

    // Set string
	result.append("fk5;Ring(");
	result.append(gammalib::str(m_centre.ra_deg()));
	result.append(",");
	result.append(gammalib::str(m_centre.dec_deg()));
	result.append(",");
	result.append(gammalib::str(m_radius1));
	result.append(",");
	result.append(gammalib::str(m_radius2));
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
 * @brief Print ring region
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing region information.
 ***************************************************************************/
std::string GSkyRegionRing::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append string
    	result.append("=== GSkyRegionRing ===");
    	result.append("\n(");
        result.append(gammalib::str(m_centre.ra_deg()));
        result.append(",");
        result.append(gammalib::str(m_centre.dec_deg()));
        result.append(",");
        result.append(gammalib::str(m_radius1));
        result.append(",");
        result.append(gammalib::str(m_radius2));
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
 * centre is not larger than the region ring (between inner and outer radius).
 ***************************************************************************/
bool GSkyRegionRing::contains(const GSkyDir& dir) const
{
	// Initialise return value
	bool dir_is_in = false;

	// calculate distance from dir to centre
	double distance = dir.dist_deg(m_centre);

	// Check if distance < radius
	if (distance <= radius2() && distance >= radius1()) {
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
bool GSkyRegionRing::contains(const GSkyRegion& reg) const
{
	// Initialise return value
	bool fully_inside = false;

	// If other region is Ring use a simple way to calculate
	if (reg.type() == "Ring") {

		// Create circular region from reg
		const GSkyRegionRing* regcirc =
              dynamic_cast<const GSkyRegionRing*>(&reg);

		// Calculate angular distance between the centres
		double ang_dist = m_centre.dist_deg(regcirc->centre());

		// Check if the region is contained in this
		if ((ang_dist + regcirc->radius2()) <= m_radius2 && (ang_dist - regcirc->radius1()) >= m_radius1) {
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
bool GSkyRegionRing::overlaps(const GSkyRegion& reg) const
{
	// Initialise return value
	bool overlap = false;

	// If other region is Ring use a simple way to calculate
	if (reg.type() == "Ring") {

		// Create circular region from reg
		const GSkyRegionRing* regcirc =
              dynamic_cast<const GSkyRegionRing*>(&reg);

		// Calculate angular distance between the centres
		double ang_dist = m_centre.dist_deg(regcirc->centre());

		// Check if the distance is smaller than the sum of both outer radii
		if (ang_dist <= (m_radius2 + regcirc->radius2())) {
			overlap = true;
		}
		
	  // Check if the distance is smaller than the sum of both inner radii
		if (ang_dist <= (m_radius1 + regcirc->radius1())) {
			overlap = true;
		}

	  // Check if two regions overlap
		if (ang_dist >= (m_radius1 + regcirc->radius1()) && ang_dist <= (m_radius2 + regcirc->radius2())) {
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
void GSkyRegionRing::init_members(void)
{
    // Initialise members
	m_centre  = GSkyDir();
	m_radius1 = 0.0;
	m_radius2 = 0.0;
	m_type    = "Ring";

	//Return
	return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] region Circular sky region.
 ***************************************************************************/
void GSkyRegionRing::copy_members(const GSkyRegionRing& region)
{
	// Copy attributes
	m_centre  = region.m_centre;
	m_radius1 = region.m_radius1;
	m_radius2 = region.m_radius2;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GSkyRegionRing::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute solid angle
 ***************************************************************************/
void GSkyRegionRing::compute_solid_angle(void)
{
	// Convert inner and outer radius from degrees to radians
	double radius1_rad = m_radius1 * gammalib::deg2rad;
	double radius2_rad = m_radius2 * gammalib::deg2rad;

	// Compute solid angle
	m_solid1 = gammalib::twopi * (1.0 - std::cos(radius1_rad));
	m_solid2 = gammalib::twopi * (1.0 - std::cos(radius2_rad));

    // Return
    return;
}
