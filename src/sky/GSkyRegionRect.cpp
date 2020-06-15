/***************************************************************************
 *               GSkyRegionRect.cpp - Circular sky region class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2016 by Michael Mayer                               *
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
 * @file GSkyRegionRect.hpp
 * @brief Circular sky region implementation
 * @author Michael Mayer
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GSkyRegionRect.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_WIDTH                            "GSkyRegionRect::width(double&)"
#define G_HEIGHT                          "GSkyRegionRect::height(double&)"
#define G_READ                         "GSkyRegionRect::read(std::string&)"
#define G_CONTAINS                  "GSkyRegionRect::contains(GSkyRegion&)"
#define G_OVERLAPS                  "GSkyRegionRect::overlaps(GSkyRegion&)"

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
GSkyRegionRect::GSkyRegionRect(void) : GSkyRegion()
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
 * @param[in] w Region width [deg].
 * @param[in] h Region height [deg].
 ***************************************************************************/
GSkyRegionRect::GSkyRegionRect(const GSkyDir& centre, const double& w,
                                   const double& h) :
                  GSkyRegion()
{
    // Initialise members
	init_members();

	// Set members
	this->centre(centre);
    this->width(w);
    this->height(h);

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
GSkyRegionRect::GSkyRegionRect(const double& ra, const double& dec,
                                   const double& w, const double& h) :
                  GSkyRegion()
{
    // Initialise members
	init_members();

	// Set members
	this->centre(ra, dec);
    this->width(w);
    this->height(h);

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
GSkyRegionRect::GSkyRegionRect(const std::string& line) : GSkyRegion()
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
 * @param[in] region Circular sky region.
 ***************************************************************************/
GSkyRegionRect::GSkyRegionRect(const GSkyRegionRect& region) :
                  GSkyRegion(region)
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
GSkyRegionRect::~GSkyRegionRect(void)
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
 * @param[in] region Circular sky region.
 * @return Circular sky region.
 ***************************************************************************/
GSkyRegionRect& GSkyRegionRect::operator=(const GSkyRegionRect& region)
{
    // Execute only if object is not identical
    if (this != &region) {

        // Copy base class members
        this->GSkyRegion::operator=(region);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(region);

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
void GSkyRegionRect::clear(void)
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
 * @brief Clone rectangular sky region
 *
 * @return Deep copy to rectangular sky region.
 ***************************************************************************/
GSkyRegionRect* GSkyRegionRect::clone(void) const
{
    // Clone rectangular sky region
    return new GSkyRegionRect(*this);
}


/***********************************************************************//**
 * @brief Set width of rectangular region
 *
 * @param[in] width Width [deg].
 *
 * @exception GException::invalid_argument
 *            Width value is less than 0.
 *
 * Sets the width of the rectangular sky region. Only non-negative radii are
 * allowed.
 ***************************************************************************/
void GSkyRegionRect::width(const double& width)
{
	if (width < 0.0) {
        std::string msg =
            "A width of "+gammalib::str(width)+" degrees has been"
            " specified for a rectangular sky region.\n"
            "The width of a rectangular sky region can't be less than 0"
            " degrees.";
		throw GException::invalid_argument(G_WIDTH, msg);
	}

    // Set radius value
    m_halfwidth = 0.5*width;

    // Compute the solid angle
    compute_solid_angle();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set height of rectangular region
 *
 * @param[in] height Height [deg].
 *
 * @exception GException::invalid_argument
 *            Height value is less than 0.
 *
 * Sets the height of the rectangular sky region. Only non-negative radii are
 * allowed.
 ***************************************************************************/
void GSkyRegionRect::height(const double& height)
{
    if (height < 0.0) {
        std::string msg =
            "A height of "+gammalib::str(height)+" degrees has been"
            " specified for a rectangular sky region.\n"
            "The height of a rectangular sky region can't be less than 0"
            " degrees.";
        throw GException::invalid_argument(G_HEIGHT, msg);
    }

    // Set radius value
    m_halfheight = 0.5*height;

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
void GSkyRegionRect::read(const std::string& line)
{
	// Clear the current instance
	clear();

	// Split the string into 2 parts seperated by #
	std::vector<std::string> substrings = gammalib::split(line,"#");
	std::string region_def = (!substrings.empty())   ? substrings[0] : "";
	std::string comment    = (substrings.size() > 1) ? substrings[1] : "";

	// Finding the box
	if (region_def.find("box") == std::string::npos) {
        std::string msg =
            "Unable to find the key word \"box\" in provided string"
            " \""+line+"\".\n"
            "The \"box\" key word is mandatory.";
		throw GException::invalid_value(G_READ, msg);
	}

	// Get the the coordinate system of the values
	std::string system = gammalib::split(region_def, ";")[0];

	// Get the substring of the important values
	unsigned    pos          = region_def.find("box(");
	unsigned    end          = region_def.find(")");
	std::string circlestring = region_def.substr(pos+7, end);
	circlestring.erase(circlestring.find(")"), 1);

	// Get the values of the region x,y,and radius
	std::vector<std::string> values = gammalib::split(circlestring,",");
	if (values.size() != 5) {
        std::string msg =
            "Invalid number of "+gammalib::str(values.size())+" arguments"
            " after the \"circle\" key word in provided string \""+line+"\".\n"
            "Exactly 5 arguments are expected.";
		throw GException::invalid_value(G_READ, msg);
	}
	double x      = gammalib::todouble(values[0]);
	double y      = gammalib::todouble(values[1]);
	double w      = gammalib::todouble(values[2]);
    double h      = gammalib::todouble(values[3]);
    double posang = gammalib::todouble(values[4]);

    // Assert posang is zero
    if (std::abs(posang)>0.1) {
        std::string msg =
            "Box position angle is not zero! Handling of rotated rectangular"
            " sky regions is generally possible but not implemented.";
        throw GException::invalid_value(G_READ, msg);
    }

    // Get radius units
    if (gammalib::contains(values[2], "'")) {
        w /= 60.0;
    }
    else if (gammalib::contains(values[2], "\"")) {
        w /= 3600.0;
    }
    if (gammalib::contains(values[3], "'")) {
        h /= 60.0;
    }
    else if (gammalib::contains(values[3], "\"")) {
        h /= 3600.0;
    }

	// Initialise centre direction
	GSkyDir centre = GSkyDir();

	// Set the values in the correct system
	if (system == "fk5" || system == "icrs") {
		centre.radec_deg(x,y);
	}
	else if (system == "galactic") {
		centre.lb_deg(x,y);
	}
	else {
        std::string msg =
            "Unsupported coordinate system \""+system+"\" in provided string"
            " \""+line+"\".\n"
            "Only the following coordinate systems are supported: "
            "\"fk5\", \"icrs\" and \"galactic\".";
		throw GException::invalid_value(G_READ, msg);
	}

	// Set members
	this->centre(centre);
    this->width(w);
    this->height(h);

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
std::string GSkyRegionRect::write(void) const
{
    // Allocate string
	std::string result;

    // Set string
	result.append("fk5;box(");
	result.append(gammalib::str(m_centre.ra_deg()));
	result.append(",");
	result.append(gammalib::str(m_centre.dec_deg()));
	result.append(",");
	result.append(gammalib::str(2*m_halfwidth));
    result.append(",");
    result.append(gammalib::str(2*m_halfheight));
    result.append(",");
    result.append(gammalib::str(0));
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
 * @brief Print rectangular region
 *
 * @param[in] chatter Chattiness
 * @return String containing region information.
 ***************************************************************************/
std::string GSkyRegionRect::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
    	result.append("=== GSkyRegionRect ===");

        // Append sky circle information
        result.append("\n"+gammalib::parformat("Right Ascension of centre"));
        result.append(gammalib::str(m_centre.ra_deg())+" deg");
        result.append("\n"+gammalib::parformat("Declination of centre"));
        result.append(gammalib::str(m_centre.dec_deg())+" deg");
        result.append("\n"+gammalib::parformat("Width"));
        result.append(gammalib::str(2*m_halfwidth)+" deg");
        result.append("\n"+gammalib::parformat("Height"));
        result.append(gammalib::str(2*m_halfheight)+" deg");

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
bool GSkyRegionRect::contains(const GSkyDir& dir) const
{
	// Initialise return value
	bool dir_is_in = false;

    // std::cout<<"GSkyRegionRect::contains: (1/2) dec: "<<m_centre.dec_deg()-m_halfheight<<"<="<<dir.dec_deg()<<"<="<<m_centre.dec_deg()+m_halfheight<<std::endl;

    // Check dec axis
    if ((dir.dec_deg()<=(m_centre.dec_deg() + m_halfheight)) &&
        (dir.dec_deg()>=(m_centre.dec_deg() - m_halfheight))) {

        // Compute width scaled for declination of queried skydir
        double scale = std::cos(dir.dec());
        double half_width_scaled = m_halfwidth/scale;

        // Check ra axis
        if ((dir.ra_deg()<=(m_centre.ra_deg() + half_width_scaled)) &&
            (dir.ra_deg()>=(m_centre.ra_deg() - half_width_scaled))) {
            dir_is_in = true;
        }

        // std::cout<<"GSkyRegionRect::contains: (2/2)  ra: "<<m_centre.ra_deg()-half_width_scaled<<"<="<<dir.ra_deg()<<"<="<<m_centre.ra_deg()+half_width_scaled<<std::endl;
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
bool GSkyRegionRect::contains(const GSkyRegion& reg) const
{
    // Initialise return value
    bool fully_inside = false;

    throw GException::feature_not_implemented(G_CONTAINS,
              "Cannot compare Rect region type yet");

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
bool GSkyRegionRect::overlaps(const GSkyRegion& reg) const
{

	// Initialise return value
	bool overlap = false;

	throw GException::feature_not_implemented(G_OVERLAPS,
          "Cannot compare Rect region type yet");

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
void GSkyRegionRect::init_members(void)
{
    // Initialise members
	m_centre = GSkyDir();
	m_halfwidth  = 0.0;
    m_halfheight = 0.0;
	m_type   = "Rect";

	//Return
	return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] region Circular sky region.
 ***************************************************************************/
void GSkyRegionRect::copy_members(const GSkyRegionRect& region)
{
	// Copy attributes
	m_centre = region.m_centre;
	m_halfwidth  = region.m_halfwidth;
    m_halfheight = region.m_halfheight;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GSkyRegionRect::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute solid angle
 *
 * As width and height are opening angles the rectangular boundaries are
 * not aligned to spherical coordinate lines, the box keeps it opening angles.
 * Hence the solidangle will vary.
 ***************************************************************************/
void GSkyRegionRect::compute_solid_angle(void)
{
    // Compute width scaling factor from box center
    double scale = std::cos(m_centre.dec());

    // Compute solid angle
    m_solid = (2*m_halfwidth)/scale * (2*m_halfheight);

    // Return
    return;
}
