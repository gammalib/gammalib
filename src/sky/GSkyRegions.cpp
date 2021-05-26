/***************************************************************************
 *                 GSkyRegions.cpp - Sky region container class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2021 by Pierrick Martin                             *
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
* @file GSkyRegions.cpp
* @brief Sky regions container class definition
* @author Pierrick Martin
*/

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#include <fstream>
#endif
#include "GBase.hpp"
#include "GTools.hpp"
#include "GSkyDir.hpp"
#include "GSkyRegion.hpp"
#include "GSkyRegionCircle.hpp"
#include "GSkyRegionRectangle.hpp"
#include "GSkyRegions.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ACCESS                       "GSkyRegions::operator[](std::string&)"
#define G_AT                                           "GSkyRegions::at(int&)"
#define G_SET1                           "GSkyRegions::set(int&, GSkyRegion&)"
#define G_SET2                   "GSkyRegions::set(std::string&, GSkyRegion&)"
#define G_APPEND                            "GSkyRegions::append(GSkyRegion&)"
#define G_INSERT1                     "GSkyRegions::insert(int&, GSkyRegion&)"
#define G_INSERT2             "GSkyRegions::insert(std::string&, GSkyRegion&)"
#define G_REMOVE1                                  "GSkyRegions::remove(int&)"
#define G_REMOVE2                          "GSkyRegions::remove(std::string&)"
#define G_EXTEND                           "GSkyRegions::extend(GSkyRegions&)"
#define G_LOAD                                 "GSkyRegions::load(GFilename&)"
#define G_SAVE                                 "GSkyRegions::save(GFilename&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GSkyRegions::GSkyRegions(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] regions region container.
 ***************************************************************************/
GSkyRegions::GSkyRegions(const GSkyRegions& regions)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(regions);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load constructor
 *
 * @param[in] filename DS9 region file.
 *
 * Constructs region container from a DS9 region file.
 ***************************************************************************/
GSkyRegions::GSkyRegions(const GFilename& filename)
{
    // Initialise members
    init_members();

    // Load XML file
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GSkyRegions::~GSkyRegions(void)
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
 * @param[in] regions region container.
 * @return region container.
 ***************************************************************************/
GSkyRegions& GSkyRegions::operator=(const GSkyRegions& regions)
{
    // Execute only if object is not identical
    if (this != &regions) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(regions);

    } // endif: object was not identical

    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                            Public methods                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear object
 *
 * Removes all regions from the container.
 ***************************************************************************/
void GSkyRegions::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 *
 * @return Pointer to deep copy of region container
 *
 * Makes a deep copy of the region container instance.
 ***************************************************************************/
GSkyRegions* GSkyRegions::clone(void) const
{
    return new GSkyRegions(*this);
}


/***********************************************************************//**
 * @brief Return pointer to region
 *
 * @param[in] index Sky region index [0,...,size()[.
 *
 * @exception GException::out_of_range
 *            Sky region index is out of range.
 *
 * Returns a pointer to the region with the specified @p index.
 ***************************************************************************/
GSkyRegion* GSkyRegions::at(const int& index)
{
    // Throw exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Sky region index",
                                       index, size());
    }

    // Return pointer
    return m_regions[index];
}


/***********************************************************************//**
 * @brief Return pointer to region (const version)
 *
 * @param[in] index Sky region index [0,...,size()[.
 *
 * @exception GException::out_of_range
 *            Sky region index is out of range.
 *
 * Returns a const pointer to the region with the specified @p index.
 ***************************************************************************/
const GSkyRegion* GSkyRegions::at(const int& index) const
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Sky Region index",
                                       index, size());
    }

    // Return pointer
    return m_regions[index];
}


/***********************************************************************//**
 * @brief Set region in container
 *
 * @param[in] index Sky region index [0,...,size()[.
 * @param[in] region Sky region.
 * @return Pointer to deep copy of sky region.
 *
 * @exception GException::out_of_range
 *            Sky region index is out of range.
 * @exception GException::invalid_value
 *            Name of region exists already in container.
 *
 * Set sky region in the container. A deep copy of the region will be made.
 ***************************************************************************/
GSkyRegion* GSkyRegions::set(const int& index, const GSkyRegion& region)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_SET1, "Sky Region index",
                                       index, size());
    }
    #endif

    // Free existing region only if it differs from current region. This
    // prevents unintential deallocation of the argument
    if ((m_regions[index] != NULL) && (m_regions[index] != &region)) {
        delete m_regions[index];
    }

    // Assign new region by cloning
    m_regions[index] = region.clone();

    // Return pointer to region
    return m_regions[index];
}


/***********************************************************************//**
 * @brief Append region to container
 *
 * @param[in] region region.
 * @return Pointer to deep copy of region.
 *
 * @exception GException::invalid_value
 *            Name of region exists already in container.
 *
 * Appends region to the container by making a deep copy of the region and
 * storing its pointer.
 ***************************************************************************/
GSkyRegion* GSkyRegions::append(const GSkyRegion& region)
{
    // Create deep copy of region
    GSkyRegion* ptr = region.clone();

    // Append deep copy of region
    m_regions.push_back(ptr);

    // Return pointer to region
    return ptr;
}


/***********************************************************************//**
 * @brief Insert region into container
 *
 * @param[in] index Sky region index [0,...,size()[.
 * @param[in] region Sky region.
 * @return Pointer to deep copy of sky region.
 *
 * @exception GException::out_of_range
 *            Sky region index is out of range.
 * @exception GException::invalid_value
 *            Name of region exists already in container.
 *
 * Inserts a @p region into the container before the region with the specified
 * @p index.
 ***************************************************************************/
GSkyRegion* GSkyRegions::insert(const int& index, const GSkyRegion& region)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (is_empty()) {
        if (index > 0) {
            throw GException::out_of_range(G_INSERT1, "Sky Region index",
                                           index, size());
        }
    }
    else {
        if (index < 0 || index >= size()) {
            throw GException::out_of_range(G_INSERT1, "Sky Region index",
                                           index, size());
        }
    }
    #endif

    // Create deep copy of region
    GSkyRegion* ptr = region.clone();

    // Inserts deep copy of region
    m_regions.insert(m_regions.begin()+index, ptr);

    // Return pointer to region
    return ptr;
}


/***********************************************************************//**
 * @brief Remove region from container
 *
 * @param[in] index Sky region index [0,...,size()[.
 *
 * @exception GException::out_of_range
 *            Sky region index is out of range.
 *
 * Remove sky region of specified @p index from container.
 ***************************************************************************/
void GSkyRegions::remove(const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_REMOVE1, "Sky Region index",
                                       index, size());
    }
    #endif

    // Delete region
    delete m_regions[index];

    // Erase region component from container
    m_regions.erase(m_regions.begin() + index);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append region container
 *
 * @param[in] regions region container.
 *
 * Append region container to the container.
 ***************************************************************************/
void GSkyRegions::extend(const GSkyRegions& regions)
{
    // Do nothing if region container is empty
    if (!regions.is_empty()) {

        // Get size. Note that we extract the size first to avoid an
        // endless loop that arises when a container is appended to
        // itself.
        int num = regions.size();

        // Reserve enough space
        reserve(size() + num);

        // Loop over all region components and append pointers to deep copies
        for (int i = 0; i < num; ++i) {

            // Append region to container
            m_regions.push_back(regions[i]->clone());

        } // endfor: looped over all regions

    } // endif: region container was not empty

    // Return
    return;
}


/***********************************************************************//**
 * @brief Check if direction is contained in one of the regions
 *
 * @param[in] dir Sky direction.
 * @return True if the sky direction is container in one of the regions,
 *         falls otherwise.
 *
 * Checks if a sky direction is contained in one of the regions of the
 * container.
 ***************************************************************************/
bool GSkyRegions::contains(const GSkyDir& dir) const
{
    // Initialise return value
    bool overlap = false;

    // Loop over regions
    for (int i = 0; i < size(); ++i) {
        overlap = m_regions[i]->contains(dir);
        if (overlap) {
            break;
        }
    }

    // Return result
    return overlap;
}


/***********************************************************************//**
 * @brief Check if region overlaps one of the regions
 *
 * @param[in] region Sky region.
 * @return True if @p region overlaps with one of the regions in the region
 *         container, false otherwise
 *
 * Tells if region overlaps one of the regions
 ***************************************************************************/
bool GSkyRegions::overlaps(const GSkyRegion& region) const
{
    // Initialise return value
    bool overlap = false;

    // Loop over regions
    for (int i = 0; i < size(); ++i) {
        overlap = m_regions[i]->overlaps(region);
        if (overlap) {
            break;
        }
    }

    // Return result
    return overlap;
}


/***********************************************************************//**
 * @brief Check if any of the regions in two containers overlap
 *
 * @param[in] regions Sky region container.
 * @return True if one of the regions in the @p regions container overlaps
 *         with one of the regions in the container, false otherwise
 *
 * Tells if all regions in  overlaps one of the regions. Note, this method
 * returns true if ANY of the regions in the two containers overlap with
 * each other
 ***************************************************************************/
bool GSkyRegions::overlaps(const GSkyRegions& regions) const
{
    // Initialize return value
    bool overlap = false;

    // Loop over each region in the container
    for (int i = 0; i < size(); ++i) {
        overlap = regions.overlaps(*m_regions[i]);
        if (overlap) {
            break;
        }
    }

    // Return result
    return overlap;
}


/***********************************************************************//**
 * @brief Load regions from DS9 region file
 *
 * @param[in] filename DS9 region filename.
 *
 * @exception GException::file_open_error
 *            File could not be opened.
 *
 * Loads all regions from a DS9 region file.
 ***************************************************************************/
void GSkyRegions::load(const GFilename& filename)
{
    // Clear any existing regions
    clear();

    // Open file. Throw an exception if opening failed.
    std::ifstream ds9file;
    ds9file.open(filename.url().c_str());
    if (ds9file.is_open()) {

        // Loop over file lines
        std::string fileline = "";
        std::string coordsys = "galactic";
        while (ds9file.good()) {

            // Read one line
            getline(ds9file,fileline);

            // If line is a comment then continue
            if (fileline[0] == '#') {
                continue;
            }

            // Check for global definition of coordinate system
            if (std::string::npos != fileline.find("fk5")) {
                coordsys = "fk5";
            }
            else if (std::string::npos != fileline.find("icrs")) {
                coordsys = "icrs";
            }

            // If region is a circle then read circle
            if (std::string::npos != fileline.find("circle")) {

                // Create instance of GSkyRegionCircle object
                GSkyRegionCircle region;

                // If coordinate system and region defined on the same line
                if ((std::string::npos != fileline.find("fk5")) ||
                    (std::string::npos != fileline.find("icrs")) ||
                    (std::string::npos != fileline.find("galactic"))) {
                    region.read(fileline);
                    append(region);
                }

                // else, prepend the coordinate system
                else {
				    std::string newfileline = coordsys;
					newfileline.append("; ");
					newfileline.append(fileline);
					region.read(newfileline);
					append(region);
				}

			} // endif: region was circle

            // ... otherwise if region is a box then read rectangle
            else if (std::string::npos != fileline.find("box")) {

                // Create instance of GSkyRegionRectangle object
                GSkyRegionRectangle region;

                // If coordinate system and region defined on the same line
                if ((std::string::npos != fileline.find("fk5")) ||
                    (std::string::npos != fileline.find("icrs")) ||
                    (std::string::npos != fileline.find("galactic"))) {
                    region.read(fileline);
                    append(region);
                }

                // else, prepend the coordinate system
                else {
                    std::string newfileline = coordsys;
                    newfileline.append("; ");
                    newfileline.append(fileline);
                    region.read(newfileline);
                    append(region);
                }

            } // endif: region was box

        } // endwhile: looped over file

        // Close file
        ds9file.close();

        // Store filename
        m_filename = filename;

    }

    // File could not be opened
    else {
        throw GException::file_open_error(G_LOAD, filename.url());
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save regions into DS9 region file
 *
 * @param[in] filename DS9 region filename.
 *
 * @exception GException::file_open_error
 *            File could not be opened.
 *
 * Saves all regions in the container into a DS9 region file.
 ***************************************************************************/
void GSkyRegions::save(const GFilename& filename) const
{
    // Open file
    std::ofstream ds9file;
    ds9file.open(filename.url().c_str());

    // If file opened correctly, then save regions
    if (ds9file.is_open()) {

        // Write global definition
        std::string fileline;
        fileline.append("# Region file format: DS9 version 4.1\n");
        fileline.append("global color=green dashlist=8 3 width=1");
        ds9file << fileline << "\n";

        // Loop over regions in container
        for (int i = 0; i < size(); ++i) {
            ds9file << m_regions[i]->write() << "\n";
        }

        // Close file
        ds9file.close();

        // Store filename
        m_filename = filename;

    }

    // ... otherwise, if file could not be opened then throw an exception
    else {
        throw GException::file_open_error(G_SAVE, filename.url());
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print regions
 *
 * @param[in] chatter Chattiness.
 * @return String containing region container information.
 *
 * Prints all regions into a string.
 ***************************************************************************/
std::string GSkyRegions::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GSkyRegions ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of regions"));
        result.append(gammalib::str(size()));

        // Append regions
        for (int i = 0; i < size(); ++i) {
            result.append("\n"+m_regions[i]->print(chatter));
        }

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                              Private methods                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GSkyRegions::init_members(void)
{
    // Initialise members
    m_filename.clear();
    m_regions.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] regions region container.
 *
 * Makes a copy of all class members. All regions are deep copied.
 ***************************************************************************/
void GSkyRegions::copy_members(const GSkyRegions& regions)
{
    // Copy members
    m_filename = regions.m_filename;

    // Copy regions
    for (int i = 0; i < regions.m_regions.size(); ++i) {
        m_regions.push_back((regions.m_regions[i]->clone()));
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 *
 * Deallocates all regions. The method loops over the region container and
 * deallocates the memory that has been allocated before.
 ***************************************************************************/
void GSkyRegions::free_members(void)
{
    // Free regions
    for (int i = 0; i < m_regions.size(); ++i) {
        if (m_regions[i] != NULL) {
            delete m_regions[i];
        }
        m_regions[i] = NULL;
    }

    // Return
    return;
}
