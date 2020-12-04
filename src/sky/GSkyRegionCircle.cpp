/***************************************************************************
 *               GSkyRegionCircle.cpp - Circular sky region class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2020 by Michael Mayer                               *
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
#include "GTools.hpp"
#include "GSkyRegionCircle.hpp"
#include "GSkyRegionRect.hpp"
#include "GSkyRegionMap.hpp"

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
GSkyRegionCircle::GSkyRegionCircle(void) : GSkyRegion()
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
 * @param[in] radius Region radius (degrees).
 ***************************************************************************/
GSkyRegionCircle::GSkyRegionCircle(const GSkyDir& centre,
                                   const double& radius) : GSkyRegion()
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
 * @param[in] ra Right Ascension of region centre (degrees).
 * @param[in] dec Declination of region centre (degrees).
 * @param[in] radius Region radius (degrees).
 ***************************************************************************/
GSkyRegionCircle::GSkyRegionCircle(const double& ra,
                                   const double& dec,
                                   const double& radius) : GSkyRegion()
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
GSkyRegionCircle::GSkyRegionCircle(const std::string& line) : GSkyRegion()
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
GSkyRegionCircle::GSkyRegionCircle(const GSkyRegionCircle& region) :
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
 * @param[in] region Circular sky region.
 * @return Circular sky region.
 ***************************************************************************/
GSkyRegionCircle& GSkyRegionCircle::operator=(const GSkyRegionCircle& region)
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
 * @param[in] radius Radius (degrees).
 *
 * @exception GException::invalid_argument
 *            Radius value is less than 0.
 *
 * Sets the radius of the circular sky region. Only non-negative radii are
 * allowed.
 ***************************************************************************/
void GSkyRegionCircle::radius(const double& radius)
{
    // Check if radius is valid
    if (radius < 0.0) {
        std::string msg =
            "A negative radius of "+gammalib::str(radius)+" degrees has been "
            "specified for a circular sky region. Please specify a "
            "non-negative radius.";
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
    std::string region_def = (!substrings.empty())   ? substrings[0] : "";
    std::string comment    = (substrings.size() > 1) ? substrings[1] : "";

    // Finding the circle
    if (region_def.find("circle") == std::string::npos) {
        std::string msg =
            "Unable to find the key word \"circle\" in provided string"
            " \""+line+"\". The \"circle\" key word is mandatory.";
        throw GException::invalid_value(G_READ, msg);
    }

    // Get the the coordinate system of the values
    std::string system = gammalib::split(region_def, ";")[0];

    // Get the substring of the important values
    size_t      pos          = region_def.find("circle(");
    size_t      end          = region_def.find(")");
    std::string circlestring = region_def.substr(pos+7, end);
    circlestring.erase(circlestring.find(")"), 1);

    // Get the values of the region x,y,and radius
    std::vector<std::string> values = gammalib::split(circlestring,",");
    if (values.size() != 3) {
        std::string msg =
            "Invalid number of "+gammalib::str(values.size())+" arguments "
            "after the \"circle\" key word in provided string \""+line+"\". "
            "Exactly 3 arguments are expected.";
        throw GException::invalid_value(G_READ, msg);
    }
    double x      = gammalib::todouble(values[0]);
    double y      = gammalib::todouble(values[1]);
    double radius = gammalib::todouble(values[2]);

    // Get radius units
    if (gammalib::contains(values[2], "'")) {
        radius /= 60.0;
    }
    else if (gammalib::contains(values[2], "\"")) {
        radius /= 3600.0;
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
            "Unsupported coordinate system \""+system+"\" in provided string "
            "\""+line+"\". Only the following coordinate systems are "
            "supported: \"fk5\", \"icrs\" and \"galactic\".";
        throw GException::invalid_value(G_READ, msg);
    }

    // Set members
    this->centre(centre);
    this->radius(radius);

    // Compute solid angle
    compute_solid_angle();

    // Check if there is a given name for the region and set it. Strings may
    // be quoted with " or ' or {}
    pos = comment.find("text");
    if (pos != std::string::npos) {
        pos = comment.find("=", pos+4);
        if (pos != std::string::npos) {

            // Set possible quotes
            const std::string quotes_open  = "{'\"";
            const std::string quotes_close = "}'\"";

            // Set quotes found flag
            bool found = false;

            // Search for possible quotes and extract name
            for (int i = 0; i < quotes_open.size(); ++i) {
                std::string quote_open  = quotes_open.substr(i,1);
                std::string quote_close = quotes_close.substr(i,1);
                pos = comment.find(quote_open, pos+1);
                if (pos != std::string::npos) {
                    end = comment.find(quote_close, pos+1);
                    if (end != std::string::npos) {
                        int length = end-pos-1;
                        if (length > 0) {
                            m_name = comment.substr(pos+1, length);
                        }
                        else {
                            m_name = "";
                        }
                        found = true;
                    }
                    else {
                        std::string msg = "Missing "+quote_open+" following "+
                                          quote_close+" encountered after "
                                          "text attribute. Please add closing "
                                          "quotes to text="+quote_open+"Name"+
                                          quote_close+" attribute.";
                        throw GException::invalid_value(G_READ, msg);
                    }
                }
            } // endfor: looped over possible quotes

            // Throw an exception if no quotes were found
            if (!found) {
                std::string msg = "Missing quotes \" or ' or {} following "
                                  "text attribute. Please specify quotes "
                                  "around the names of the text, for example "
                                  "text={Name}.";
                throw GException::invalid_value(G_READ, msg);
            }

        } // endif: = character found

        // Signal that no = character was found after text
        else {
            std::string msg = "Missing = following text attribute.  Please "
                              "specify the value for a text attribute, "
                              "for example text={Name}.";
            throw GException::invalid_value(G_READ, msg);
        }

    } // endif: text attribute found

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
        result.append("{");
        result.append(m_name);
        result.append("}");
    }

    // Return string
    return result;
}


/***********************************************************************//**
 * @brief Print circular region
 *
 * @param[in] chatter Chattiness
 * @return String containing region information.
 ***************************************************************************/
std::string GSkyRegionCircle::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GSkyRegionCircle ===");

        // Append sky circle information
        result.append("\n"+gammalib::parformat("Right Ascension of centre"));
        result.append(gammalib::str(ra())+" deg");
        result.append("\n"+gammalib::parformat("Declination of centre"));
        result.append(gammalib::str(dec())+" deg");
        result.append("\n"+gammalib::parformat("Radius"));
        result.append(gammalib::str(radius())+" deg");

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
    bool contains = false;

    // If the other region is a circle then check whether it is fully
    // container in the circle
    if (reg.type() == "Circle") {

        // Create circular region from reg
        const GSkyRegionCircle* circle = dynamic_cast<const GSkyRegionCircle*>(&reg);

        // Calculate angular distance between the centres
        double ang_dist = m_centre.dist_deg(circle->centre());

        // Check if the region is contained in this
        if ((ang_dist + circle->radius()) <= m_radius) {
            contains = true;
        }

    } // endif: other region was of type "Circle"

    // ... otherwise, if the other region is a rectangle then check whether
    // all four corners are contained within the circle
    else if (reg.type() == "Rect") {

        // Create rectangular region from reg
        const GSkyRegionRect* rect = dynamic_cast<const GSkyRegionRect*>(&reg);

        // Check containment of all corners
        contains = this->contains(rect->corner(0)) &&
                   this->contains(rect->corner(1)) &&
                   this->contains(rect->corner(2)) &&
                   this->contains(rect->corner(3));

    }

    // ... otherwise, if the other region is a map then check whether it
    // is contained within the rectangle
    else if (reg.type() == "Map") {

        // Create map from reg
        const GSkyRegionMap* map = dynamic_cast<const GSkyRegionMap*>(&reg);

        // Get non-zero indices
        std::vector<int> indices = map->nonzero_indices();

        // Initialise containment flag to true
        contains = true;

        // Loop over all indices
        for (int i = 0; i < indices.size(); ++i) {

            // Get sky direction of map pixel
            GSkyDir dir = map->map().inx2dir(i);

            // If pixel is not contained then set containment flag to false
            // and break
            if (!this->contains(dir)) {
                contains = false;
                break;
            }

        } // endfor: looped over non-zero map indices

    } // endif: other region was of type "Map"

    // ... otherwise throw an exception
    else {
        throw GException::feature_not_implemented(G_CONTAINS,
              "Cannot compare to region type \""+reg.type()+"\" yet.");
    }

    // Return containment flag
    return contains;
}


/***********************************************************************//**
 * @brief Checks if region is overlapping with this region
 *
 * @param[in] reg Sky region.
 *
 * @exception GException::feature_not_implemented
 *            Regions differ in type.
 *
 * @todo Implement checks for rectangles and maps
 ***************************************************************************/
bool GSkyRegionCircle::overlaps(const GSkyRegion& reg) const
{
    // Initialise return value
    bool overlap = false;

    // If region is circle then check overlap by comparing the centre
    // distances
    if (reg.type() == "Circle") {

        // Create circular region from reg
        const GSkyRegionCircle* circle = dynamic_cast<const GSkyRegionCircle*>(&reg);

        // Calculate angular distance between the centres
        double ang_dist = m_centre.dist_deg(circle->centre());

        // Check if the distance is smaller than the sum of both radii
        if (ang_dist <= (m_radius + circle->radius())) {
            overlap = true;
        }

    } // endif: region was circle

    // ... otherwise if region is rectangle then check overlap
    else if (reg.type() == "Rect") {

        // Create rectangular region from reg
        const GSkyRegionRect* rect = dynamic_cast<const GSkyRegionRect*>(&reg);

        // Check overlap with circle
        overlap = rect->overlaps(*this);

    } // endif: region was rectangle

    // ... otherwise if region is map then check overlap with map
    else if (reg.type() == "Map") {

        // Create map from reg
        const GSkyRegionMap* map = dynamic_cast<const GSkyRegionMap*>(&reg);

        // Check overlap with circle
        overlap = map->overlaps(*this);

    } // endif: region was map

    // ... otherwise throw an exception
    else {
        throw GException::feature_not_implemented(G_CONTAINS,
              "Cannot compare to region type \""+reg.type()+"\" yet.");
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


/***********************************************************************//**
 * @brief Equality operator
 *
 * @param[in] a First sky region circle.
 * @param[in] b Second sky region circle.
 * @return True if both sky region circles are identical.
 *
 * Returns true if two sky region circles and identical.
 ***************************************************************************/
bool operator==(const GSkyRegionCircle &a, const GSkyRegionCircle &b)
{
    // Return equality
    return ((a.m_centre == b.m_centre) && (a.m_radius == b.m_radius));
}


/***********************************************************************//**
 * @brief Non equality operator
 *
 * @param[in] a First sky region circle.
 * @param[in] b Second sky region circle.
 * @return True if both sky region circles are different.
 ***************************************************************************/
bool operator!=(const GSkyRegionCircle &a, const GSkyRegionCircle &b)
{
    // Return non equality
    return ((a.m_centre != b.m_centre) || (a.m_radius != b.m_radius));
}
