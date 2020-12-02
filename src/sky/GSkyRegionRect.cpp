/***************************************************************************
 *              GSkyRegionRect.cpp - Rectangular sky region class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2019-2020 by Andreas Specovius                           *
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
 * @brief Rectangular sky region implementation
 * @author Andreas Specovius
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//#include "GMath.hpp"
#include "GTools.hpp"
#include "GSkyRegionCircle.hpp"
#include "GSkyRegionRect.hpp"
#include "GSkyRegionMap.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_WIDTH                              "GSkyRegionRect::width(double&)"
#define G_HEIGHT                            "GSkyRegionRect::height(double&)"
#define G_READ                           "GSkyRegionRect::read(std::string&)"
#define G_CONTAINS                    "GSkyRegionRect::contains(GSkyRegion&)"
#define G_OVERLAPS                    "GSkyRegionRect::overlaps(GSkyRegion&)"
#define G_GET_CORNER                       "GSkyRegionRect::get_corner(int&)"

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
 * @brief Sky direction constructor
 *
 * @param[in] centre Centre of rectangle.
 * @param[in] width Region width (degrees).
 * @param[in] height Region height (degrees).
 * @param[in] posang Position angle (degrees).
 ***************************************************************************/
GSkyRegionRect::GSkyRegionRect(const GSkyDir& centre,
                               const double&  width,
                               const double&  height,
                               const double&  posang) : GSkyRegion()
{
    // Initialise members
    init_members();

    // Set members
    this->centre(centre);
    this->width(width);
    this->height(height);
    this->posang(posang);

    // Compute solid angle
    compute_solid_angle();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Right Ascension and Declination constructor
 *
 * @param[in] ra Right Ascension of region centre [deg].
 * @param[in] dec Declination of region centre [deg].
 * @param[in] width Region width (degrees).
 * @param[in] height Region height (degrees).
 * @param[in] posang Position angle (degrees).
 ***************************************************************************/
GSkyRegionRect::GSkyRegionRect(const double& ra,
                               const double& dec,
                               const double& width,
                               const double& height,
                               const double& posang) : GSkyRegion()
{
    // Initialise members
    init_members();

    // Set members
    this->centre(ra, dec);
    this->width(width);
    this->height(height);
    this->posang(posang);

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
 * @param[in] region Rectangular sky region.
 ***************************************************************************/
GSkyRegionRect::GSkyRegionRect(const GSkyRegionRect& region) : GSkyRegion(region)
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
 * @param[in] region Rectangular sky region.
 * @return Rectangular sky region.
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
 * @param[in] width Width (degrees).
 *
 * @exception GException::invalid_argument
 *            Width value is less than 0.
 *
 * Sets the @p width of the rectangular sky region. Only non-negative widths
 * are allowed.
 ***************************************************************************/
void GSkyRegionRect::width(const double& width)
{
    // Throw an exception if the width is less than zero
    if (width < 0.0) {
        std::string msg =
            "A negative width of "+gammalib::str(width)+" degrees has been "
            "specified for a rectangular sky region. Please specify a "
            "non-negative width.";
        throw GException::invalid_argument(G_WIDTH, msg);
    }

    // Set radius value
    m_width = width;

    // Compute the solid angle
    compute_solid_angle();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set height of rectangular region
 *
 * @param[in] height Height (degrees).
 *
 * @exception GException::invalid_argument
 *            Height value is less than 0.
 *
 * Sets the @p height of the rectangular sky region. Only non-negative
 * heights are allowed.
 *
 * Note that for a position angle of zero, the height axis is pointing to
 * celestial north.
 ***************************************************************************/
void GSkyRegionRect::height(const double& height)
{
    if (height < 0.0) {
        std::string msg =
            "A negative height of "+gammalib::str(height)+" degrees has been "
            "specified for a rectangular sky region. Please specify a "
            "non-negative height.";
        throw GException::invalid_argument(G_HEIGHT, msg);
    }

    // Set radius value
    m_height = height;

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
            "Unable to find the key word \"box\" in provided string "
            "\""+line+"\". The \"box\" key word is mandatory.";
        throw GException::invalid_value(G_READ, msg);
    }

    // Get the the coordinate system of the values
    std::string system = gammalib::split(region_def, ";")[0];

    // Get the substring of the important values
    unsigned    pos          = region_def.find("box(");
    unsigned    end          = region_def.find(")");
    std::string regionstring = region_def.substr(pos+4, end);
    regionstring.erase(regionstring.find(")"), 1);

    // Get the values of the region x,y,and radius
    std::vector<std::string> values = gammalib::split(regionstring,",");
    if (values.size() != 5) {
        std::string msg =
            "Invalid number of "+gammalib::str(values.size())+" arguments "
            "after the \"box\" key word in provided string \""+line+"\". "
            "Exactly 5 arguments are expected.";
        throw GException::invalid_value(G_READ, msg);
    }
    double x      = gammalib::todouble(values[0]);
    double y      = gammalib::todouble(values[1]);
    double width  = gammalib::todouble(values[2]);
    double height = gammalib::todouble(values[3]);
    double posang = gammalib::todouble(values[4]);

    // Get radius units
    if (gammalib::contains(values[2], "'")) {
        width /= 60.0;
    }
    else if (gammalib::contains(values[2], "\"")) {
        width /= 3600.0;
    }
    if (gammalib::contains(values[3], "'")) {
        height /= 60.0;
    }
    else if (gammalib::contains(values[3], "\"")) {
        height /= 3600.0;
    }
    if (gammalib::contains(values[4], "'")) {
        posang /= 60.0;
    }
    else if (gammalib::contains(values[4], "\"")) {
        posang /= 3600.0;
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
    this->width(width);
    this->height(height);
    this->posang(posang);

    // Compute solid angle
    compute_solid_angle();

    // Check if there is a given name for the region and set it
    std::vector<std::string>comments = gammalib::split(comment, " ");
    for (int i = 0; i < comments.size(); ++i) {
        if (gammalib::contains(comments[i], "text")) {
            std::vector<std::string> attributes = gammalib::split(comments[i], "=");
            if (attributes.size() < 2) {
                std::string msg =
                    "Invalid character sequence encountered in provided "
                    "string \""+line+"\". An attribute of the type "
                    "\"text=Name\" is expected.";
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
    result.append(gammalib::str(m_width));
    result.append(",");
    result.append(gammalib::str(m_height));
    result.append(",");
    result.append(gammalib::str(m_posang));
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
 * @brief Checks if sky direction lies within region
 *
 * @param[in] dir Sky direction.
 *
 * A sky direction lies within a region when its distance to the region
 * centre is not larger than the region extension in both axes directions.
 ***************************************************************************/
bool GSkyRegionRect::contains(const GSkyDir& dir) const
{
    // Compute sky direction in local coordinate system
    GSkyDir locdir = transform_to_local(dir);

    // Check containment using local coordinate system
    bool dir_is_in = contains_local(locdir);

    // Return bool
    return dir_is_in;
}


/***********************************************************************//**
 * @brief Checks if local direction lies within region
 *
 * @param[in] locdir Sky direction in local coordinate system.
 *
 * A local direction lies within a region when its distance to the region
 * centre is not larger than the region extension in both axes directions.
 ***************************************************************************/
bool GSkyRegionRect::contains_local(const GSkyDir& locdir) const
{
    // Initialise return value
    bool dir_is_in = false;

    // Check containment: declination axis
    if ((std::abs(locdir.dec_deg()) - 1.0e-10) <= 0.5 * m_height) {

        // Check containment: right ascension axis
        if ((std::abs(locdir.ra_deg()) - 1.0e-10) <= 0.5 * m_width) {

            // Sky direction is inside the rectangle
            dir_is_in = true;
        }
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
 *            Not all region types supported currently.
 ***************************************************************************/
bool GSkyRegionRect::contains(const GSkyRegion& reg) const
{
    // Initialise return value
    bool is_fully_inside = false;

    // If other region is circle use a simple way to calculate
    if (reg.type() == "Circle") {

        // Create circular region from reg
        const GSkyRegionCircle* regcirc =
              dynamic_cast<const GSkyRegionCircle*>(&reg);

        // Transform circle center to local coordinate system
        GSkyDir local_centre = transform_to_local(regcirc->centre());

        // Get circle radius
        const double circ_rad = regcirc->radius();

        // Check extension of circle along declination axis
        if ((std::abs(local_centre.dec_deg()) + circ_rad) <= 0.5 * m_height) {

            // Check extension of circle along right ascension axis
            if ((std::abs(local_centre.ra_deg()) + circ_rad) <= 0.5 * m_width) {

                // Circle is fully contained in this rectangle
                is_fully_inside = true;
            }
        }
    } // Region was of type "Circle"

    // If other region is rectangle use a simple way to calculate
    else if (reg.type() == "Rect") {

        // Create rectangular region from reg
        const GSkyRegionRect* regrect =
              dynamic_cast<const GSkyRegionRect*>(&reg);

        // Loop over the four corners of the rectangle :regrect:
        for (int icorner = 0; icorner < 4; ++icorner) {

            // Get skydir of current corner
            GSkyDir corner = regrect->get_corner(icorner);

            // Check corner containment
            if (!contains(corner)) {
                break;
            }
            // Did we arrive at the last corner?
            else if (icorner == 3) {

                // As we reached this point, the three prior corners were inside
                // and the last corner now is also inside:
                // Rectangle is fully inside
                is_fully_inside = true;
            }
        } // Looped over the four corners
    } // Region was of type "Rect"

    // ... otherwise throw an exception
    else {
        throw GException::feature_not_implemented(G_CONTAINS,
              "Cannot compare rectangular region with other than rectangle or"
              " circle yet.");
    }

    // Return value
    return is_fully_inside;
}


/***********************************************************************//**
 * @brief Checks if region is overlapping with this region
 *
 * @param[in] reg Sky region.
 *
 * @exception GException::feature_not_implemented
 *            Regions differ in type.
 *
 * @todo: - Improve implementation for rectangle-rectangle
 ***************************************************************************/
bool GSkyRegionRect::overlaps(const GSkyRegion& reg) const
{
    // Initialise return value
    bool is_overlapping = false;

    // If other region is circle use a simple way to calculate
    if (reg.type() == "Circle") {

        // Create circular region from reg
        const GSkyRegionCircle* regcirc =
              dynamic_cast<const GSkyRegionCircle*>(&reg);

        // Transform circle center to local coordinate system
        GSkyDir local_centre = transform_to_local(regcirc->centre());

        // Get circle radius
        const double circ_rad = regcirc->radius();

        // Check extension of circle along declination axis
        if (std::abs(local_centre.dec_deg()) <= (0.5 * m_height + circ_rad)) {

            // Check extension of circle along right ascension axis
            if (std::abs(local_centre.ra_deg()) <= (0.5 * m_width + circ_rad)) {

                // Circle is fully contained in this rectangle
                is_overlapping = true;
            }
        }
    } // Region was of type "Circle"

    else if (reg.type() == "Rect") {

        // Create rectangular region from reg
        const GSkyRegionRect* regrect =
              dynamic_cast<const GSkyRegionRect*>(&reg);

        // Dirty kludge: compare vs map
        GSkyRegionMap regmap = GSkyRegionMap(regrect);
        is_overlapping = regmap.overlaps(*this);

    } // Region was of type "Rect"

    // ... otherwise throw an exception
    else {
        throw GException::feature_not_implemented(G_OVERLAPS,
              "Cannot compare rectangular region with other than rectangle or"
              " circle yet.");
    }

    // Return value
    return is_overlapping;
}


/***********************************************************************//**
 * @brief Transform a sky direction to the local cartesian coordinate system.
 *
 * @param[in] skydir Sky direction in global coordinate system.
 * @return Sky direction in local cartesian region object coordinates.
 *
 * Transform the sky direction @p skydir to the local cartesian coordinate
 * system. The origin of the local coordinate system is fixed to the center
 * of the rectangle and aligned in +ra (width) and +dec (height). Please
 * note, that a GSkyDir object is returned, even if the local coordsys is
 * is cartesian.
 ***************************************************************************/
GSkyDir GSkyRegionRect::transform_to_local(const GSkyDir& skydir) const
{
    // Get distance and polar angle relative to rectangle centre
    double dist = m_centre.dist(skydir);
    double pa   = m_centre.posang(skydir) - m_posang * gammalib::deg2rad;

    // Compute corresponding x and y coordinate from polar coords
    double x = dist * std::sin(pa);
    double y = dist * std::cos(pa);

    // Create output sky direction
    GSkyDir transformed;
    transformed.radec(x, y);

    // Return
    return transformed;
}


/***********************************************************************//**
 * @brief Transform a local coordinate to the global coordinate system.
 *
 * @param[in] locdir Direction in local coordinate system.
 * @return Sky direction in global coordinate system.
 *
 * Transform the local direction @p locdir to the global coordinate system.
 * The origin of the local coordinate system is fixed to the center of the
 * rectangle and aligned in +ra (width) and +dec (height).
 ***************************************************************************/
GSkyDir GSkyRegionRect::transform_to_global(const GSkyDir& locdir) const
{
    // // Initialise centre in local coordinates (maybe make class member?)
    GSkyDir centre_local;
    centre_local.radec(0,0);

    // Access x and y in local coordsys
    double x = locdir.ra();
    double y = locdir.dec();

    // Compute dist and position angle in polar coords
    double dist = std::sqrt(x*x + y*y);
    double pa   = gammalib::pihalf - std::atan2(y, x) +
                  m_posang  * gammalib::deg2rad;

    // Create output global sky direction
    GSkyDir transformed(m_centre);
    transformed.rotate(pa, dist);

    // Return
    return transformed;
}


/***********************************************************************//**
 * @brief Get the sky direction of a corner of the rectangle, there are four.
 *
 * @param[in] index Corner index [0,3]
 * @return Sky direction of corner in global coordinate system.
 *
 * @exception GException::out_of_range
 *            Corner index is not in [0,3].
 *
 * The index counts CCW starting from +RA,+DEC: +R+D, +R-D, -R-D, -R+D.
 ***************************************************************************/
GSkyDir GSkyRegionRect::get_corner(const int& index) const
{
    // Assert index is in [0,3]
    if ((index < 0) || (index > 3)) {
        throw GException::out_of_range(G_GET_CORNER, "Corner index", index, 3);
    }

    // Compute the offset to the corners
    // Index counts CCW starting from +RA,+DEC: +R+D, +R-D, -R-D, -R+D
    double dx = 0.5 * m_width  * (1 - 2*int(index >= 2));
    double dy = 0.5 * m_height * (1 - 2*int((index==1) || (index==2)));

    // Define local corner sky direction
    GSkyDir corner;
    corner.radec_deg(dx,dy);

    // Transform to global sky direction
    corner = transform_to_global(corner);

    // Return
    return corner;
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

        // Append sky region information
        result.append("\n"+gammalib::parformat("Right Ascension of centre"));
        result.append(gammalib::str(ra())+" deg");
        result.append("\n"+gammalib::parformat("Declination of centre"));
        result.append(gammalib::str(dec())+" deg");
        result.append("\n"+gammalib::parformat("Width"));
        result.append(gammalib::str(width())+" deg");
        result.append("\n"+gammalib::parformat("Height"));
        result.append(gammalib::str(height())+" deg");
        result.append("\n"+gammalib::parformat("PA"));
        result.append(gammalib::str(posang())+" deg");

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
void GSkyRegionRect::init_members(void)
{
    // Initialise members
    m_type   = "Rect";
    m_centre.clear();
    m_width  = 0.0;
    m_height = 0.0;
    m_posang = 0.0;

    //Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] region Rectangular sky region.
 ***************************************************************************/
void GSkyRegionRect::copy_members(const GSkyRegionRect& region)
{
    // Copy attributes
    m_centre = region.m_centre;
    m_width  = region.m_width;
    m_height = region.m_height;
    m_posang = region.m_posang;

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
 * @brief Compute solid angle (sr)
 ***************************************************************************/
void GSkyRegionRect::compute_solid_angle(void)
{
    // Compute solid angle
    m_solid = m_width * m_height * gammalib::deg2rad * gammalib::deg2rad;

    // Return
    return;
}
