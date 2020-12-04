/***************************************************************************
 *                GSkyRegionMap.cpp - Sky region map class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2020 by Pierrick Martin                             *
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
 * @file GSkyRegionMap.cpp
 * @brief Sky region map implementation
 * @author Pierrick Martin
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GMath.hpp"
#include "GSkyDir.hpp"
#include "GSkyMap.hpp"
#include "GSkyRegion.hpp"
#include "GSkyRegionCircle.hpp"
#include "GSkyRegionRectangle.hpp"
#include "GSkyRegionMap.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_CONTAINS                     "GSkyRegionMap::contains(GSkyRegion&)"

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
GSkyRegionMap::GSkyRegionMap(void) : GSkyRegion()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief FITS filename constructor
 *
 * @param[in] filename FITS file.
 *
 * Constructs region from a FITS file.
 ***************************************************************************/
GSkyRegionMap::GSkyRegionMap(const GFilename& filename) : GSkyRegion()
{
    // Initialise members
    init_members();

    // Load FITS file
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Sky map constructor
 *
 * @param[in] map Sky map.
 ***************************************************************************/
GSkyRegionMap::GSkyRegionMap(const GSkyMap& map) : GSkyRegion()
{
    // Initialise members
    init_members();

    // Set sky map
    this->map(map);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Region conversion constructor
 *
 * @param[in] region Pointer to sky region.
 *
 * Constructs a sky region map from any region type.
 ***************************************************************************/
GSkyRegionMap::GSkyRegionMap(const GSkyRegion* region) : GSkyRegion(*region)
{
    // Initialise members
    init_members();

    // Get dynamic casts
    const GSkyRegionMap*       map    = dynamic_cast<const GSkyRegionMap*>(region);
    const GSkyRegionCircle*    circle = dynamic_cast<const GSkyRegionCircle*>(region);
    const GSkyRegionRectangle* rect   = dynamic_cast<const GSkyRegionRectangle*>(region);

    // If region is of type GSkyRegionMap then copy class members
    if (map != NULL) {
        copy_members(*map);
    }

    // ... otherwise, if region is of type GSkyRegionCircle then set region
    // map from region circle
    else if (circle != NULL) {
        set_region_circle(circle);
    }

    // If region is of type GSkyRegionRectangle than set region map fron
    // region rect
    else if (rect != NULL) {
        set_region_rectangle(rect);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] region Sky region map.
 ***************************************************************************/
GSkyRegionMap::GSkyRegionMap(const GSkyRegionMap& region) : GSkyRegion(region)
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
GSkyRegionMap::~GSkyRegionMap(void)
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
 * @param[in] region Sky region map.
 * @return Sky region map.
 ***************************************************************************/
GSkyRegionMap& GSkyRegionMap::operator=(const GSkyRegionMap& region)
{
    // Execute only if object is not identical
    if (this != &region) {

        // Free members
        free_members();

        // Initialise private members
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
void GSkyRegionMap::clear(void)
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
 * @brief Clone sky region map
 *
 * @return Deep copy to sky region map.
 ***************************************************************************/
GSkyRegionMap* GSkyRegionMap::clone(void) const
{
    // Clone sky map region
    return new GSkyRegionMap(*this);
}


/***********************************************************************//**
 * @brief Set sky map
 *
 * @param[in] map Sky map
 *
 * Sets the sky map of a region map.
 ***************************************************************************/
void GSkyRegionMap::map(const GSkyMap& map)
{
    // Set map object
    m_map = map;

    // Set non-zero pixel indices
    set_nonzero_indices();

    // Compute solid angle
    compute_solid_angle();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load region map from FITS file
 *
 * @param[in] filename FITS file name.
 ***************************************************************************/
void GSkyRegionMap::load(const GFilename& filename)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Load map
    m_map.load(filename);

    // Set non-zero indices
    set_nonzero_indices();

    // Compute solid angle
    compute_solid_angle();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Create region from a string in DS9 format
 *
 * @param[in] line String describing region in DS9 region format.
 *
 * Empty implementation of a virtual function because map is not supported
 * as DS9 region.
 *
 * @todo Translate any DS9 region into a map.
 ***************************************************************************/
void GSkyRegionMap::read(const std::string& line)
{
    // Do nothing because map is not supported as DS9 region
    return;
}


/***********************************************************************//**
 * @brief Write string describing region in DS9 region format
 *
 * @return DS9 region string
 *
 * Returns an empty string.
 *
 * @todo Translate region map into a DS9 polygon.
 ***************************************************************************/
std::string GSkyRegionMap::write(void) const
{
    // Allocate string
    std::string result;

    // Return string
    return result;
}


/***********************************************************************//**
 * @brief Print region description
 *
 * @param[in] chatter Chattiness.
 * @return String containing region information.
 ***************************************************************************/
std::string GSkyRegionMap::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append string
        result.append("=== GSkyRegionMap ===");
        result.append("\n(");
        result.append(m_map.print());
        result.append(")");

    } // endif: chatter was not silent

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Check if a given direction is contained in this region
 *
 * @param[in] dir Sky direction.
 * @return True if sky direction is in region, false otherwise.
 *
 * Checks if sky direction is contained in the region map.
 ***************************************************************************/
bool GSkyRegionMap::contains(const GSkyDir& dir) const
{
    // Initialize
    bool contain = false;

    // If direction is within map boundaries
    if (m_map.contains(dir)) {

        // Convert sky direction into pixel index
        int i = m_map.dir2inx(dir);

        // Set containment flag
        contain = (m_map(i,0) != 0.0);

    } // endif: map contains direction

    // Return value
    return contain;
}


/***********************************************************************//**
 * @brief Checks if a given region is fully contained within this region
 *
 * @param[in] reg Sky region.
 * @return True if the region is contained within this region, false
 *         otherwise.
 *
 * @exception GException::feature_not_implemented
 *            Specified sky region type not supported by method.
 ***************************************************************************/
bool GSkyRegionMap::contains(const GSkyRegion& reg) const
{
    // Initialise return value
    bool inside = true;

    // Case of a region map
    if (reg.type() == "Map") {

        // Cast to map type
        const GSkyRegionMap* inreg = static_cast<const GSkyRegionMap*>(&reg);

        // Retrieve map data
        const GSkyMap&          regmap   = inreg->map();
        const std::vector<int>& regnzvec = inreg->nonzero_indices();

        // Loop over input map non-zero pixels and check if they are all in
        for (int i = 0; i < regnzvec.size(); ++i) {
            if (!contains(regmap.inx2dir(regnzvec[i]))) {
                inside = false;
                break;
            }
        }

    } // endcase: region map

    // Case of a circular region
    else if (reg.type() == "Circle") {

        // Create circular region from reg
        const GSkyRegionCircle* inreg =
              static_cast<const GSkyRegionCircle*>(&reg);

        // Retrieve circle data
        const GSkyDir& centre = inreg->centre();
        double         rad    = inreg->radius();

        // Test a number of points along the circle
        const int numtest  = 360;
        double    rotangle = 360.0 / numtest;
        GSkyDir   testdir;
        for (int i = 0; i < numtest; ++i) {
            testdir = centre;
            testdir.rotate_deg(i*rotangle, rad);
            if (!contains(testdir)) {
                inside = false;
                break;
            }
        }

    } // endcase: circular region

    // Case of a rectangular region
    else if (reg.type() == "Rectangle") {

        // Create rectangular region from reg
        const GSkyRegionRectangle* rect =
              static_cast<const GSkyRegionRectangle*>(&reg);

        // Create region map from rectangle
        GSkyRegionMap rect_map(rect);

        // Retrieve map data
        const GSkyMap&          regmap   = rect_map.map();
        const std::vector<int>& regnzvec = rect_map.nonzero_indices();

        // Loop over input map non-zero pixels and check if they are all in
        for (int i = 0; i < regnzvec.size(); ++i) {
            if (!contains(regmap.inx2dir(regnzvec[i]))) {
                inside = false;
                break;
            }
        }
    } // endcase: rectangular region

    // Other cases not implemented
    else {
        throw GException::feature_not_implemented(G_CONTAINS,
            "Method only implemented for cicular, rectangular and map regions.");
    }

    // Return result
    return inside;
}


/***********************************************************************//**
 * @brief Checks if a given region is overlapping with this region
 *
 * @param[in] reg Sky region.
 * @return True if the region is overlapping, false otherwise.
 ***************************************************************************/
bool GSkyRegionMap::overlaps(const GSkyRegion& reg) const
{
    // Initialise return value
    bool overlap = false;

    // Loop over non-zero pixels until direction is contained into region
    for (int i = 0; i < m_nonzero_indices.size(); ++i) {
        GSkyDir dir = m_map.inx2dir(m_nonzero_indices[i]);
        if (reg.contains(dir)) {
            overlap = true;
            break;
        }
    }

    // Return
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
void GSkyRegionMap::init_members(void)
{
    // Initialie base class members
    m_type = "Map";

    // Initialise members
    m_map.clear();
    m_nonzero_indices.clear();

    //Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] map Sky region map.
 ***************************************************************************/
void GSkyRegionMap::copy_members(const GSkyRegionMap& map)
{
    // Copy members
    m_map             = map.m_map;
    m_nonzero_indices = map.m_nonzero_indices;

    // Compute solid angle
    compute_solid_angle();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GSkyRegionMap::free_members(void)
{
    // No memory to be freed, return
    return;
}


/***********************************************************************//**
 * @brief Compute solid angle
 ***************************************************************************/
void GSkyRegionMap::compute_solid_angle(void)
{
    // Initialise solid angle
    m_solid = 0.0;

    // Loop over all non-zero pixels of map and add up solid angles
    for (int i = 0; i < m_nonzero_indices.size(); ++i) {
        m_solid += m_map.solidangle(m_nonzero_indices[i]);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Create an array of non-zero pixel indices
 ***************************************************************************/
void GSkyRegionMap::set_nonzero_indices(void)
{
    // Clear zero pixel array
    m_nonzero_indices.clear();

    // Loop over map pixels and find non-zero pixels
    for (int i = 0; i < m_map.npix(); ++i) {
        if (m_map(i,0) != 0.0) {
            m_nonzero_indices.push_back(i);
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set region map from region circle
 *
 * @param[in] circle Region circle
 *
 * Sets a region map from a region circle. The method allocates a sky map
 * with a default resolution of 1/20 of the circle radius. In any case, the
 * sky map resolution is comprised within [0.0002, 0.1] degrees. The default
 * size of the sky map is 2.5 times the circle radius, and in any case, at
 * least 10 sky map pixels will be allocated. The sky map will be in
 * celestial coordinates and TAN projection.
 ***************************************************************************/
void GSkyRegionMap::set_region_circle(const GSkyRegionCircle* circle)
{
    // Set sky map resolution from circle radius and constrain the
    // resolution to [0.0002, 0.1] deg
    double resolution = 0.05 * circle->radius();
    if (resolution < 0.0002) {
        resolution = 0.0002;
    }
    else if (resolution > 0.1) {
        resolution = 0.1;
    }

    // Set number of sky map pixels
    int nbins = int(2.5 * circle->radius() / resolution);
    if (nbins < 10) {
        nbins = 10;
    }

    // Set sky map
    GSkyMap map("TAN", "CEL", circle->ra(), circle->dec(),
                              resolution,   resolution,
                              nbins,        nbins);

    // Set sky map pixels
    for (int i = 0; i < map.npix(); ++i) {

        // Compute distance to circle centre
        double distance = circle->centre().dist_deg(map.inx2dir(i));

        // If distance is smaller or equal to radius then set map
        // pixel to 1
        if (distance <= circle->radius()) {
            map(i) = 1.0;
        }

    } // endfor: looped over sky map pixels

    // Set sky map as region map
    this->map(map);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set region map from region rectangle
 *
 * @param[in] rect Region rectangle
 *
 * Sets a region map from a region rectangle.
 *
 * A skymap is allocated with a resolution goal of 0.01 degrees (same for
 * both axes). The particular resolution of either axes depends on the width
 * over height of the RoI of the rotated rectangle. The map extent is chosen
 * dynamically to match the RoI of the rotated rectangle. Three sparse pixels
 * are added at either side of the map. The sky map is created in celestial
 * coordinates and TAN projection.
 ***************************************************************************/
void GSkyRegionMap::set_region_rectangle(const GSkyRegionRectangle* rect)
{
    // Define desired resolution (in degrees)
    const double resolution_goal = 0.01;

    // Try to estimate extent of the (rotated) rectangle:
    // Loop over the rectangle corners to get the extension
    double dxmax = 0.0;
    double dymax = 0.0;
    for (int icorner = 0; icorner < 4; ++icorner) {

        // Get skydir of current corner
        GSkyDir corner = rect->corner(icorner);

        // Compute polar props in global coords relative to rectangle centre
        double dist = rect->centre().dist(corner);
        double pa   = rect->centre().posang(corner);

        // Compute required map extension along RA/Dec for this corner
        double dx = dist * std::sin(pa);
        double dy = dist * std::cos(pa);

        // Select the largest extension
        if (dy > dymax) {
            dymax = dy;
        }
        if (dx > dxmax) {
            dxmax = dx;
        }

    } // endfor: looped over all corners

    // Compute width of map RoI in degrees
    double map_width  = 2.0 * dxmax * gammalib::rad2deg;
    double map_height = 2.0 * dymax * gammalib::rad2deg;

    // Compute num of bins in x direction, round
    int nbins_x = int(map_width  / resolution_goal + 0.5);
    int nbins_y = int(map_height / resolution_goal + 0.5);

    // Ensure a minimum num of bins
    if (nbins_x < 10) {
        nbins_x = 10;
    }
    if (nbins_y < 10) {
        nbins_y = 10;
    }

    // Compute actual resolution using nbins
    double resolution_x = map_width  / nbins_x;
    double resolution_y = map_height / nbins_y;

    // Add leakage pixels
    nbins_x += 6;
    nbins_y += 6;

    // Create fine sky map
    GSkyMap map_fine("TAN", "CEL", rect->ra(),   rect->dec(),
                                  -resolution_x, resolution_y,
                                   nbins_x,      nbins_y);

    // Set sky map pixels that are contained in rectangle to 1
    for (int ipix = 0; ipix < map_fine.npix(); ++ipix) {
        if (rect->contains(map_fine.inx2dir(ipix))) {
            map_fine(ipix) = 1.0;
        }
    }

    // Set sky map as region map
    this->map(map_fine);

    // Return
    return;
}
