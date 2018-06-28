/***************************************************************************
 *                  GSkyRegions.hpp - Sky region container class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2017 by Pierrick Martin                             *
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
 * @file GSkyRegions.hpp
 * @brief Sky regions container class definition
 * @author Pierrick Martin
 */

#ifndef GSKYREGIONS_HPP
#define GSKYREGIONS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GContainer.hpp"
#include "GSkyRegion.hpp"
#include "GException.hpp"
#include "GFilename.hpp"

/* __ Forward declarations _______________________________________________ */
class GSkyDir;


/***********************************************************************//**
 * @class GSkyRegions
 *
 * @brief Sky region container class
 *
 * This container class collects sky regions objects derived from the
 * GSkyRegion abstract class. These can be accessed from the access 
 * operator[] or through the at() method with index checking. Access to 
 * the region can be done from index or name, so the name of each region 
 * in the container has to be unique.
 *
 * The object can be initialised from a DS9 region file, and has load and 
 * save methods from/to a DS9 region file.
 ***************************************************************************/
class GSkyRegions : public GContainer {

public:
    // Constructors and destructors
    GSkyRegions(void);
    GSkyRegions(const GSkyRegions& regions);
    explicit GSkyRegions(const GFilename& filename);
    virtual ~GSkyRegions(void);

    // Operators
    GSkyRegions&      operator=(const GSkyRegions& regions);
    GSkyRegion*       operator[](const int& index);
    const GSkyRegion* operator[](const int& index) const;

    // Methods
    void              clear(void);
    GSkyRegions*      clone(void) const;
    std::string       classname(void) const;
    GSkyRegion*       at(const int& index);
    const GSkyRegion* at(const int& index) const;
    int               size(void) const;
    bool              is_empty(void) const;
    GSkyRegion*       set(const int& index, const GSkyRegion& region);
    GSkyRegion*       append(const GSkyRegion& region);
    GSkyRegion*       insert(const int& index, const GSkyRegion& region);
    void              remove(const int& index);
    void              reserve(const int& num);
    void              extend(const GSkyRegions& regions);
    bool              contains(const GSkyDir& dir) const;
    bool              overlaps(const GSkyRegion& region) const;
    bool              overlaps(const GSkyRegions& regions) const;
    void              load(const GFilename& filename);
    void              save(const GFilename& filename) const;
    const GFilename&  filename(void) const;
    void              filename(const GFilename& filename);
    std::string       print(const GChatter& chatter = NORMAL) const;
	
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GSkyRegions& regions);
    void free_members(void);

    // Protected members
    mutable GFilename        m_filename;   //!< Filename of origin
    std::vector<GSkyRegion*> m_regions;    //!< List of regions
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GSkyRegions").
 ***************************************************************************/
inline
std::string GSkyRegions::classname(void) const
{
    return ("GSkyRegions");
}


/***********************************************************************//**
 * @brief Return pointer to region
 *
 * @param[in] index region index [0,...,size()-1].
 *
 * Returns a pointer to the region with the specified @p index.
 ***************************************************************************/
inline
GSkyRegion* GSkyRegions::operator[](const int& index)
{
    return (m_regions[index]);
}


/***********************************************************************//**
 * @brief Return pointer to region (const version)
 *
 * @param[in] index region index [0,...,size()-1].
 *
 * Returns a const pointer to the region with the specified @p index.
 ***************************************************************************/
inline
const GSkyRegion* GSkyRegions::operator[](const int& index) const
{
    return (m_regions[index]);
}


/***********************************************************************//**
 * @brief Return number of regions in container
 *
 * @return Number of regions in container.
 *
 * Returns the number of regions in the region container.
 ***************************************************************************/
inline
int GSkyRegions::size(void) const
{
    return (int)m_regions.size();
}


/***********************************************************************//**
 * @brief Signals if there are no regions in container
 *
 * @return True if container is empty, false otherwise.
 *
 * Signals if the region container does not contain any region.
 ***************************************************************************/
inline
bool GSkyRegions::is_empty(void) const
{
    return (m_regions.empty());
}


/***********************************************************************//**
 * @brief Reserves space for regions in container
 *
 * @param[in] num Number of regions
 *
 * Reserves space for @p num regions in the container.
 ***************************************************************************/
inline
void GSkyRegions::reserve(const int& num)
{
    m_regions.reserve(num);
    return;
}


/***********************************************************************//**
 * @brief Return regions file name
 *
 * @return File name from which the regions have been loaded or into which
 *         the regions have been saved.
 *
 * Returns the file name from which the regions have been loaded or into
 * which the regions have been saved. The returned string will be set if
 * the load(), save() or filename() setter methods were called.
 ***************************************************************************/
inline
const GFilename& GSkyRegions::filename(void) const
{
    return (m_filename);
}


/***********************************************************************//**
 * @brief Set regions file name
 *
 * @param[in] filename Regions file name.
 *
 * Set the regions file name.
 ***************************************************************************/
inline
void GSkyRegions::filename(const GFilename& filename)
{
    m_filename = filename;
    return;
}

#endif /* GSKYREGIONS_HPP */
