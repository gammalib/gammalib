/***************************************************************************
 *                 GSkyRegions.cpp - Sky region container class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Pierrick Martin                                  *
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
#include "GSkyRegion.hpp"
#include "GSkyRegionCircle.hpp"
#include "GSkyRegions.hpp"
#include "GTools.hpp"

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
#define G_LOAD                               "GSkyRegions::load(std::string&)"
#define G_SAVE                               "GSkyRegions::save(std::string&)"

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
GSkyRegions::GSkyRegions(const std::string& filename)
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


/***********************************************************************//**
 * @brief Return pointer to region
 *
 * @param[in] name region name.
 *
 * @exception GException::invalid_argument
 *            region with specified name not found in container.
 *
 * Returns a pointer to the region with the specified @p name.
 ***************************************************************************/
GSkyRegion* GSkyRegions::operator[](const std::string& name)
{
    // Get region index
    int index = get_index(name);

    // Throw exception if region name was not found
    if (index == -1) {
		std::string msg ="Region with name \""+name+"\" not found in container.";
        throw GException::invalid_argument(G_ACCESS, msg);
    }

    // Return pointer
    return m_regions[index];
}


/***********************************************************************//**
 * @brief Return pointer to region (const version)
 *
 * @param[in] name region name.
 *
 * @exception GException::invalid_argument
 *            region with specified name not found in container.
 *
 * Returns a const pointer to the region with the specified @p name.
 ***************************************************************************/
const GSkyRegion* GSkyRegions::operator[](const std::string& name) const
{
    // Get region index
    int index = get_index(name);

    // Throw exception if region name was not found
    if (index == -1) {
		std::string msg ="Region with name \""+name+"\" not found in container.";
        throw GException::invalid_argument(G_ACCESS, msg);		
    }

    // Return pointer
    return m_regions[index];
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
 * @param[in] index region index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            region index is out of range.
 *
 * Returns a pointer to the region with the specified @p index.
 ***************************************************************************/
GSkyRegion* GSkyRegions::at(const int& index)
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, index, 0, size()-1);
    }

    // Return pointer
    return m_regions[index];
}


/***********************************************************************//**
 * @brief Return pointer to region (const version)
 *
 * @param[in] index region index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            region index is out of range.
 *
 * Returns a const pointer to the region with the specified @p index.
 ***************************************************************************/
const GSkyRegion* GSkyRegions::at(const int& index) const
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, index, 0, size()-1);
    }

    // Return pointer
    return m_regions[index];
}


/***********************************************************************//**
 * @brief Set region in container
 *
 * @param[in] index sky region index [0,...,size()-1].
 * @param[in] region sky region.
 * @return Pointer to deep copy of sky region.
 *
 * @exception GException::out_of_range
 *            region index is out of range.
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
        throw GException::out_of_range(G_SET1, index, 0, size()-1);
    }
    #endif

    // Check if a region with specified name does not yet exist
    int inx = get_index(region.name());
    if (inx != -1 && inx != index) {
        std::string msg =
            "Attempt to set region with name \""+region.name()+"\" in region"
            " container at index "+gammalib::str(index)+", but a region with"
            " the same name exists already at index "+gammalib::str(inx)+
            " in the container.\n"
			"Every region in the region container needs a unique name.";
        throw GException::invalid_value(G_SET1, msg);		
    }

    // Delete any existing region
    if (m_regions[index] != NULL) delete m_regions[index];

    // Assign new region by cloning
    m_regions[index] = region.clone();

    // Return pointer to region
    return m_regions[index];
}


/***********************************************************************//**
 * @brief Set region in container
 *
 * @param[in] name region name.
 * @param[in] region region pointer.
 * @return Pointer to deep copy of region.
 *
 * @exception GException::invalid_argument
 *            region with specified name not found in container.
 * @exception GException::invalid_value
 *            Name of region exists already in container.
 *
 * Set region in the container. A deep copy of the region will be made.
 ***************************************************************************/
GSkyRegion* GSkyRegions::set(const std::string& name, const GSkyRegion& region)
{
    // Get parameter index
    int index = get_index(name);

    // Throw exception if parameter name was not found
    if (index == -1) {
        throw GException::invalid_argument(G_SET2, name);
    }

    // Check if a region with specified name does not yet exist
    int inx = get_index(region.name());
    if (inx != -1 && inx != index) {
        std::string msg =
            "Attempt to set region with name \""+region.name()+"\" in region"
            " container at index "+gammalib::str(index)+", but a region with"
            " the same name exists already at index "+gammalib::str(inx)+
            " in the container.\n"
            "Every region in the region container needs a unique name.";
        throw GException::invalid_value(G_SET2, msg);
    }

    // Delete any existing region
    if (m_regions[index] != NULL) delete m_regions[index];

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
    // Check if a region with specified name does not yet exist
    int inx = get_index(region.name());
    if (inx != -1) {
        std::string msg = 
            "Attempt to append region with name \""+region.name()+"\" to region"
            " container, but a region with the same name exists already at"
            " index "+gammalib::str(inx)+" in the container.\n"
            "Every region in the region container needs a unique name.";
        throw GException::invalid_value(G_APPEND, msg);
    }

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
 * @param[in] index region index [0,...,size()-1].
 * @param[in] region region.
 * @return Pointer to deep copy of region.
 *
 * @exception GException::out_of_range
 *            region index is out of range.
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
    if (isempty()) {
        if (index > 0) {
            throw GException::out_of_range(G_INSERT1, index, 0, size()-1);
        }
    }
    else {
        if (index < 0 || index >= size()) {
            throw GException::out_of_range(G_INSERT1, index, 0, size()-1);
        }
    }
    #endif

    // Check if a region with specified name does not yet exist
    int inx = get_index(region.name());
    if (inx != -1) {
        std::string msg =
            "Attempt to insert region with name \""+region.name()+"\" in region"
            " container before index "+gammalib::str(index)+", but a region"
            " with the same name exists already at index "+gammalib::str(inx)+
            " in the container.\n"
            "Every region in the region container needs a unique name.";
        throw GException::invalid_value(G_INSERT1, msg);
    }

    // Create deep copy of region
    GSkyRegion* ptr = region.clone();

    // Inserts deep copy of region
    m_regions.insert(m_regions.begin()+index, ptr);

    // Return pointer to region
    return ptr;
}


/***********************************************************************//**
 * @brief Insert region into container
 *
 * @param[in] name region name.
 * @param[in] region region.
 * @return Pointer to deep copy of region.
 *
 * @exception GException::invalid_argument
 *            region with specified name not found in container.
 * @exception GException::invalid_value
 *            Name of region exists already in container.
 *
 * Inserts a @p region into the container before the region with the specified
 * @p name.
 ***************************************************************************/
GSkyRegion* GSkyRegions::insert(const std::string& name, const GSkyRegion& region)
{
    // Get parameter index
    int index = get_index(name);

    // Throw exception if parameter name was not found
    if (index == -1) {
        throw GException::invalid_argument(G_INSERT2, name);
    }

    // Check if a region with specified name does not yet exist
    int inx = get_index(region.name());
    if (inx != -1) {
        std::string msg =
            "Attempt to insert region with name \""+region.name()+"\" in region"
            " container before index "+gammalib::str(index)+", but a region"
            " with the same name exists already at index "+gammalib::str(inx)+
            " in the container.\n"
            "Every region in the region container needs a unique name.";
        throw GException::invalid_value(G_INSERT2, msg);
    }

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
 * @param[in] index region index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            region index is out of range.
 *
 * Remove region of specified @p index from container.
 ***************************************************************************/
void GSkyRegions::remove(const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_REMOVE1, index, 0, size()-1);
    }
    #endif

    // Erase region component from container
    m_regions.erase(m_regions.begin() + index);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Remove region from container
 *
 * @param[in] name region name.
 *
 * @exception GException::invalid_argument
 *            region with specified name not found in container.
 *
 * Remove region of specified @p name from container.
 ***************************************************************************/
void GSkyRegions::remove(const std::string& name)
{
    // Get parameter index
    int index = get_index(name);

    // Throw exception if region name was not found
    if (index == -1) {
        throw GException::invalid_argument(G_REMOVE2, name);
    }

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
    if (!regions.isempty()) {

        // Get size. Note that we extract the size first to avoid an
        // endless loop that arises when a container is appended to
        // itself.
        int num = regions.size();

        // Reserve enough space
        reserve(size() + num);

        // Loop over all region components and append pointers to deep copies 
        for (int i = 0; i < num; ++i) {

            // Check if region name does not yet exist
            int inx = get_index(regions[i]->name());
            if (inx != -1) {
                std::string msg =
                    "Attempt to append region with name \""+regions[i]->name()+
                    "\" to region container, but a region with the same name"
                    " exists already at index "+gammalib::str(inx)+" in the"
                    " container.\n"
                    "Every region in the region container needs a unique name.";
                throw GException::invalid_value(G_EXTEND, msg);
            }

            // Append region to container
            m_regions.push_back(regions[i]->clone());

        } // endfor: looped over all regions

    } // endif: region container was not empty
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Signals if region name exists
 *
 * @param[in] name region name.
 * @return True if region with specified @p name exists.
 *
 * Searches all region names for a match with the specified @p name. If the
 * specified name has been found, true is returned.
 ***************************************************************************/
bool GSkyRegions::contains(const std::string& name) const
{
    // Get region index
    int index = get_index(name);

    // Return
    return (index != -1);
}


/***********************************************************************//**
 * @brief Load regions from DS9 region file
 *
 * @param[in] filename DS9 region filename.
 *
 * @exception GException::file_open_error
 *            File could not be opened.
 *
 * Loads all regions from a DS9 region file. See the read() method for more
 * information about the expected structure of the XML file.
 ***************************************************************************/
void GSkyRegions::load(const std::string& filename)
{
    // Clear any existing regions
    clear();

    // Open file. Throw an exception if opening failed.
    std::ifstream ds9file;
    ds9file.open(filename.c_str());
	if (ds9file.is_open()) {
        
		// Loop over file lines
		std::string fileline="";
		std::string coordsys="galactic";
        while ( ds9file.good() ) {
			
			// Read one line
			getline (ds9file,fileline);
			
			// If line is a comment then continue
			if (fileline[0] == '#') {
                continue;
            }
			
			// Check for global definition of coordinate system
			if (std::string::npos != fileline.find("fk5")) {
				coordsys="fk5";
			}
			
			// If region is a circle
			if (std::string::npos != fileline.find("circle")) {

				// Create instance of GSkyRegion object
				GSkyRegionCircle region;				

				// If coordinate system and region defined on the same line
				if ((std::string::npos != fileline.find("fk5")) ||
					(std::string::npos != fileline.find("galactic"))) {
					region.read(fileline);
					append(region);
				}
                
				// else, prepend the coordinate system
                else {
				    std::string newfileline=coordsys;
					newfileline.append("; ");
					newfileline.append(fileline);
					region.read(newfileline);
					append(region);
				}
			} 
		}
		
        // Close file
        ds9file.close();

        // Store filename
        m_filename = filename;

	} 

	// File could not be opened
    else {
        throw GException::file_open_error(G_LOAD, filename);
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
void GSkyRegions::save(const std::string& filename) const
{
    // Open file
    std::ofstream ds9file;
    ds9file.open(filename.c_str());

    // If file opened correctly, then save regions
    if (ds9file.is_open()) {
    
		// Write global definition
		std::string fileline;
		fileline.append("# Region file format: DS9 version 4.1\n");
		fileline.append("global color=green dashlist=8 3 width=1");
//		fileline.append("font=\"helvetica 10 normal\" select=1");
//                                "highlite=1 dash=0 fixed=0 edit=1 move=1"+
//                                "delete=1 include=1 source=1";
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
        throw GException::file_open_error(G_SAVE, filename);
    }	
		
	// Return
    return;
}


/***********************************************************************//**
 * @brief Print regions
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
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


/***********************************************************************//**
 * @brief Tells if direction is contained in one of the regions
 *
 * @param[in] dir A sky direction
 * @return True or False
 *
 * Tells if direction is contained in one of the regions
 ***************************************************************************/
 bool GSkyRegions::contains(const GSkyDir& dir) const
 {
	 // Initialise return value
	 bool dir_is_in = false;
	 
	 // Loop over regions
	 for (int i = 0; i < size(); ++i) {
		 dir_is_in=m_regions[i]->contains(dir);
		 if (dir_is_in) break;
	 }
	 
	 // Return result
	 return dir_is_in;
 }

 
/***********************************************************************//**
 * @brief Tells if region overlaps one of the regions
 *
 * @param[in] reg A sky region
 * @return True or False
 *
 * Tells if region overlaps one of the regions
 ***************************************************************************/
  bool GSkyRegions::overlaps(const GSkyRegion& reg) const
  {
	  // Initialise return value
	  bool reg_is_in = false;
	  
	  // Loop over regions
	  for (int i = 0; i < size(); ++i) {
		  reg_is_in=m_regions[i]->overlaps(reg);
		  if (reg_is_in) break;
	  }
	  
	  // Return result
	  return reg_is_in;
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
 * Makes a copy of all class members. All regions are deep copied, and the
 * linear pointer array for parameter access through the GOptimizerPars
 * base class is set.
 ***************************************************************************/
void GSkyRegions::copy_members(const GSkyRegions& regions)
{
    // Copy members
    m_filename = regions.m_filename;

    // Copy regions
    m_regions.clear();
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
        delete m_regions[i];
        m_regions[i] = NULL;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return region index by name
 *
 * @param[in] name region name.
 * @return region index (-1 if region name has not been found)
 *
 * Returns region index based on the specified @p name. If no region with the
 * specified @p name is found the method returns -1.
 ***************************************************************************/
int GSkyRegions::get_index(const std::string& name) const
{
    // Initialise index
    int index = -1;

    // Search region with specified name
    for (int i = 0; i < size(); ++i) {
        if (m_regions[i]->name() == name) {
            index = i;
            break;
        }
    }

    // Return index
    return index;
}
