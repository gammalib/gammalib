/***************************************************************************
 *               GSkyRegionMap.cpp -  Sky region map class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Anneli Schulz and Pierrick Martin                *
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
 * @author Anneli Schulz and Pierrick Martin
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GSkyRegionMap.hpp"
#include "GSkyRegionCircle.hpp"
#include "GSkyDir.hpp"
#include "GTools.hpp"

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
 * @brief Set "write" function
 *
 * @return String to be written
 ***************************************************************************/
std::string GSkyRegionMap::write(void) const
{
    // Allocate string
	std::string result;

    // Set string
	result.append("Warning: This function is not needed.");

    // Return string
	return result;
}

/***********************************************************************//**
 * @brief Direction constructor
 *
 * @param[in] centre Centre sky direction.
 * @param[in] radius Region radius [deg].
 ***************************************************************************/
GSkyRegionMap::GSkyRegionMap(const GSkyRegionMap& region)
{
    // Initialise members
	init_members();

	// Set members
	this->map(region.map());

	// Compute solid angle
	compute_solid_angle();
	
	// Get non-zero indices
	std::vector<int> pix_array;
	get_array_nonzero(pix_array);

    // Return
    return;
}

/***********************************************************************//**
 * @brief String constructor
 *
 * @param[in] filename FITS file.
 *
 * Constructs region from a FITS file.
 ***************************************************************************/
GSkyRegionMap::GSkyRegionMap(const GFits& file)
{
    // Initialise members
    init_members();

    // Load FITS file
    load(file);

    // Return
    return;
}

/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] region sky map region.
 ***************************************************************************/
GSkyRegionMap::GSkyRegionMap(const GSkymap& map)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(map);

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
 * @brief Read region from DS9 string
 *
 * @param[in] filename String in DS9 format.
 *
 * @exception GException::invalid_value
 *            Invalid value found in DS9 format string.
 ***************************************************************************/
/*void GSkyRegionMap::load(const GFits& file)
{
	// Clear the current instance
	clear();

    

	// Return
	return;
}
*/
/***********************************************************************//**
 * @brief Load skymap from FITS file.
 *
 * @param[in] filename FITS file name..
 *
 * Loads HEALPix and non HEALPix skymaps. First searches for HEALPix map in
 * FITS file by scanning all HDUs for PIXTYPE=HEALPIX. If no HEALPix map has
 * been found then search load first non-empty image.
 *
 * @todo Do we have to restrict a HEALPix map to a BinTable and a WCS map
 * to a Double precision image???
 ***************************************************************************/
void GSkyRegionMap::load(const std::string& filename)
{
    // Free memory and initialise members
    free_members();
    init_members();

    // Open FITS file
    GFits fits(filename);

    // Get number of HDUs
    int num = fits.size();

    // Initialize load flag
    bool loaded = false;

    // First search for HEALPix extension. We can skip the first extension
    // since this is always an image and a HEALPix map is stored in a
    // binary table
    for (int extno = 1; extno < num; ++extno) {

        // Get reference to HDU
        const GFitsHDU& hdu = *fits.at(extno);
        
        // If PIXTYPE keyword equals "HEALPIX" then load map
        if (hdu.has_card("PIXTYPE") && hdu.string("PIXTYPE") == "HEALPIX") {
            read_healpix(static_cast<const GFitsTable&>(hdu));
            loaded = true;
            break;
        }

    } // endfor: looped over HDUs

    // If we have not found a HEALPIX map then search now for image.
    // Skip empty images
    if (!loaded) {
        for (int extno = 0; extno < num; ++extno) {

            // Get referene to HDU
            const GFitsHDU& hdu = *fits.at(extno);

            // Skip if extension is not an image
            if (extno > 0) {
                if (hdu.string("XTENSION") != "IMAGE")
                    continue;
            }

            // Load WCS map
            read_wcs(static_cast<const GFitsImage&>(hdu));
            //loaded = true;
            break;

        } // endfor: looped over HDUs
    } // endif: no HEALPix map found

    // Close FITS file
    fits.close();

    // Return
    return;
}

/***********************************************************************//**
 * @brief Read skymap from FITS HDU
 *
 * @param[in] hdu FITS HDU.
 ***************************************************************************/
void GSkyRegionMap::read(const GFitsHDU& hdu)
{
    // Free memory and initialise members
    free_members();
    init_members();

    // Initialize load flag
    bool loaded = false;

    // If PIXTYPE keyword equals "HEALPIX" then load map
    if (hdu.has_card("PIXTYPE") && hdu.string("PIXTYPE") == "HEALPIX") {
        read_healpix(static_cast<const GFitsTable&>(hdu));
        loaded = true;
    }

    // ... otherwise try loading as non HEALPix map
    if (!loaded) {

        // Load only if HDU contains an image
        if (hdu.exttype() == 0) {
            read_wcs(static_cast<const GFitsImage&>(hdu));
            //loaded = true;
        }

    } // endif

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
/*std::string GSkyRegionMap::read(void) const
{
    // Allocate string
	std::string result;

    // Set string
	result.append("Warning: This function is not needed and virtual.");

    // Return string
	return result;
}
*/

/***********************************************************************//**
 * @brief Print circular region
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing region information.
 ***************************************************************************/
/*std::string GSkyRegionMap::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append string
    	result.append("=== GSkyRegionMap ===");
    	result.append("\n(");
        result.append(gammalib::str(m_map));
        result.append(")");

    } // endif: chatter was not silent

    // Return result
    return result;
}
*/
/***********************************************************************//**
 * @brief Tells if direction is contained in one of the region maps
 *
 * @param[in] dir A Sky direction.
 *
 * @return True or False
 *
 * Tells if direction is contained in one of the regions
 ***************************************************************************/
bool GSkyRegionMap::contains(const GSkyDir& dir) const
{
    bool fully_inside = false;
    
    // Get non zeros indices of sky region map
    std::vector<int> pix_array;
    get_array_nonzero(pix_array);
    
    // Convert sky direction into sky pixel
    GSkyPixel pixel = dir2pix(dir);
    
    for (int i = 0; i < pix_array.size(); i++ ) {
        GSkyDir direc = m_map.pix2dir(pix_array[i]);
        if (i < pix_array.size() && dir==direc) { 
            fully_inside = true;
        }
        else {
            throw GException::feature_not_implemented(G_CONTAINS,
              "Cannot compare two different region types yet");
        }
    }
    
	// Return value
	return fully_inside;
}

/***********************************************************************//**
 * @brief Checks if region is fully contained within this region
 *
 * @param[in] reg Sky region map.
 *
 * @exception GException::feature_not_implemented
 *            Regions differ in type.
 ***************************************************************************/
bool GSkyRegionMap::contains(const GSkyRegion& reg) const
{
	// Initialise return value
	bool fully_inside = false;

    // Get non zeros indices of sky region map
    std::vector<int> pix_array;
    get_array_nonzero(pix_array);
    
    for (int i = 0; i < pix_array.size(); i++ ) {
        GSkyDir dir = m_map.pix2dir(pix_array[i]);
        if (i < pix_array.size() && reg.contains(dir)) { 
            fully_inside = true;
        }
        else {
            throw GException::feature_not_implemented(G_CONTAINS,
              "Cannot compare two different region types yet");
        }
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
bool GSkyRegionMap::overlaps(const GSkyRegion& reg) const
{
	// Initialise return value
	bool overlap = false;

    // Get non zeros indices of sky region map
    std::vector<int> pix_array;
    get_array_nonzero(pix_array);
    
    for (int i = 0; i < pix_array.size(); i++ ) {
        GSkyDir dir = m_map.pix2dir(pix_array[i]);
        reg.contains(dir);
        if (overlap) break;
    }

    // Return if region is in map or not
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
	m_map  = GSkymap();
	m_type = "Map";

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

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GSkyRegionMap::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute solid angle
 ***************************************************************************/
void GSkyRegionMap::compute_solid_angle(void)
{
	// Loop over all pixels of map
	double value=0.0;
	for (int i =0; i < m_map.npix(); ++i)  {
		GSkyPixel pix=m_map.inx2pix(i);
		std::cout << pix.size() << std::endl;
		if (pix.size() > 0){
		    value += m_map.solidangle(pix);
		}
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
 * @brief Create an array of pixels
 ***************************************************************************/
void GSkyRegionMap::get_array_nonzero(std::vector<int> pix_array) const
{
    //std::vector<int> pix_array; 

    for (int i = 0; i < m_map.npix(); i++) 
    {
        GSkyPixel pix=m_map.inx2pix(i);
        std::cout << pix.size() << std::endl;
		if (pix.size() > 0){
		    pix_array.push_back(pix);
		} 
    }
    
    std::cout << pix_array.size() << std::endl;
    
}
