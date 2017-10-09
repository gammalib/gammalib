/***************************************************************************
 *               GSkyRegionMap.cpp -  Sky region map class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2016 by Pierrick Martin                             *
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
 * @file GSkyRegionMap.hpp
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
#include "GSkyRegionMap.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ                         "GSkyRegionMap::read(std::string&)"
#define G_CONTAINS                  "GSkyRegionMap::contains(GSkyRegion&)"
#define G_OVERLAPS                  "GSkyRegionMap::overlaps(GSkyRegion&)"

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
 * @brief Copy constructor
 *
 * @param[in] Map sky region
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
 * @brief Skymap constructor
 *
 * @param[in] map sky map object.
 ***************************************************************************/
GSkyRegionMap::GSkyRegionMap(const GSkyMap& map) : GSkyRegion()
{
    // Initialise members
    init_members();

    // Attach map
    m_map = map;
    
    // Get non-zero pixel indices
    get_nonzero_pixels();
    
    // Compute solid angle
    compute_solid_angle();

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
 * @param[in] sky map region.
 * @return sky map region.
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
 * @brief Clone sky map region
 *
 * @return Deep copy to sky map region.
 ***************************************************************************/
GSkyRegionMap* GSkyRegionMap::clone(void) const
{
    // Clone sky map region
    return new GSkyRegionMap(*this);
}

/***********************************************************************//**
 * @brief Load region from FITS file name
 *
 * @param[in] filename FITS file name
 ***************************************************************************/
void GSkyRegionMap::load(const GFilename& filename)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Load map
    m_map.load(filename);
    
    // Store filename
    m_filename = filename;
    
    // Get non-zero indices
    get_nonzero_pixels();
    
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
 * Empty implementation of a virtual function because.
 * map is not supported as DS9 region
 ***************************************************************************/
void GSkyRegionMap::read(const std::string& line)
{
    // Do nothing because map is not supported as DS9 region
    return;
}

/***********************************************************************//**
 * @brief Write string describing region in DS9 region format
 *
 * @return string
 ***************************************************************************/
std::string GSkyRegionMap::write(void) const
{
    // Allocate string
    std::string result;

    // Set string
    result.append("Warning: This function cannot be implemented for GSkyRegionMap.");

    // Return string
    return result;
}

/***********************************************************************//**
 * @brief Print region description
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 *
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
 * @param[in] dir A Sky direction.
 *
 * @return True or False
 *
 * Checks if direction is contained in this region
 ***************************************************************************/
bool GSkyRegionMap::contains(const GSkyDir& dir) const
{
    bool contain = false;
    
    // Convert sky direction into pixel index
    int i_dir = m_map.dir2inx(dir);
    
    for (int i = 0; i < m_nzarray.size(); i++ ) {
        if (i_dir == m_nzarray[i]) { 
            contain = true;
            break;
        }
    } // Looped over all non-zero map pixels
    
    // Return value
    return contain;
}

/***********************************************************************//**
 * @brief Checks if a given region is fully contained within this region
 *
 * @param[in] reg Sky region.
 *
 * @return True or False
 *
 * @exception GException::feature_not_implemented
 *            Regions differ in type.
 ***************************************************************************/
bool GSkyRegionMap::contains(const GSkyRegion& reg) const
{
    // Initialise return value
    bool inside = true;
    
    // Case of a map region
    if (reg.type() == "Map") {
        
		// Cast to map type
		const GSkyRegionMap* inreg =
              dynamic_cast<const GSkyRegionMap*>(&reg);
        
        // Retrieve map data
        const GSkyMap& regmap=inreg->map();
        std::vector<int> regnzvec=inreg->nonzeroindices();
        int regnznum=regnzvec.size();
        int reginnum=0;
        
        // Loop over input map non-zero pixels and check if they are all in
        for (int i = 0; i < regnznum; i++ ) {
            if (contains(regmap.inx2dir(regnzvec[i]))) {
                reginnum += 1;
            }
        }
        
        // If all non-zero pixels of map are in, region is fully contained
        if (reginnum != regnznum) {
            inside = false;
        }
    
    // Case of a circle region
    } else if (reg.type() == "Circle") {
        
		// Create circular region from reg
		const GSkyRegionCircle* inreg =
              dynamic_cast<const GSkyRegionCircle*>(&reg);
        // Retrieve circle data
        const GSkyDir& cendir=inreg->centre();
        double rad=inreg->radius();
        
        // Test a number of points along the circle
        int numtest=100;
        double rotangle=180./numtest;
        GSkyDir testdir;
        for (int i = 0; i < numtest; i++ ) {
            testdir=cendir;
            testdir.rotate_deg(i*rotangle, rad);
            if (!contains(testdir)) {
                inside = false;
                break;
            }
        }
    
    // Other cases not implemented
    } else {
            throw GException::feature_not_implemented(G_CONTAINS,
              "Only implemented for GSkyRegionMap and GSkyRegionCircle.");
    }
    
    // Return value
    return inside;
}

/***********************************************************************//**
 * @brief Checks if a given region is overlapping with this region
 *
 * @param[in] reg Sky region.
 *
 * @return True or False
 ***************************************************************************/
bool GSkyRegionMap::overlaps(const GSkyRegion& reg) const
{
    // Initialise return value
    bool overlap = false;
    
    // Loop over non-zero pixels until direction is contained into region 
    for (int i = 0; i < m_nzarray.size(); i++ ) {
        GSkyDir dir = m_map.inx2dir(m_nzarray[i]);
        if (reg.contains(dir)) {
            overlap = true;
            break;
        }
    } // Looped over non-zero map pixels

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
    // Initialise members
    m_map  = GSkyMap();
    m_type = "Map";
    m_solid = 0.0;
    m_nznum = 0;
    m_filename.clear();

    //Return
    return;
}

/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] region sky map region.
 ***************************************************************************/
void GSkyRegionMap::copy_members(const GSkyRegionMap& map)
{
    // Copy attributes
    m_map = map.m_map;
    m_filename = map.m_filename;
    
    // Get non-zero pixel indices
    get_nonzero_pixels();
    
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
    // Loop over all non-zero pixels of map and add up solid angles
    double value=0.0;
    for (int i =0; i < m_nzarray.size(); ++i) {
        value += m_map.solidangle(m_nzarray[i]);
    }
    
    // Update member value
    m_solid=0.0;
    if (!gammalib::is_infinite(value))  {
        m_solid = value;
    }
    
    // Return
    return;
}

/***********************************************************************//**
 * @brief Create an array of non-zero pixel indices
 ***************************************************************************/
void GSkyRegionMap::get_nonzero_pixels(void)
{
    // Loop over map pixels and find non-zero pixels
    for (int i = 0; i < m_map.npix(); i++) 
    {
        if (m_map(i,0) > 0) {
            m_nzarray.push_back(i);
        }
    } // Looped over all map pixels
    
    // Number of non-zero pixels
    m_nznum=m_nzarray.size();
    
    // Return
    return;
}
