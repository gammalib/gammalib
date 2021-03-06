/***************************************************************************
 *               GSkyRegionMap.hpp - Sky region map class                  *
 * ----------------------------------------------------------------------- *
 * copyright (C) 2017-2020 by Pierrick Martin                              *
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
 * @brief Sky region map class interface definition
 * @author Pierrick Martin
 */

#ifndef GSKYREGIONMAP_HPP
#define GSKYREGIONMAP_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GSkyRegion.hpp"
#include "GSkyMap.hpp"
#include "GFilename.hpp"

/* __ Forward declarations _______________________________________________ */
class GSkyDir;
class GSkyRegionMap;
class GSkyRegionRectangle;


/***********************************************************************//**
 * @class GSkyRegionMap
 *
 * @brief Interface for a sky region map
 *
 * This class provides an implementation for a sky region defined by a
 * sky map. The map is provided as a FITS file. Pixels with non-zero values
 * are within the region, pixels with zero values outside the region.
 ***************************************************************************/
class GSkyRegionMap : public GSkyRegion {

public:
    // Constructors and destructors
    GSkyRegionMap(void);
    explicit GSkyRegionMap(const GFilename& filename);
    explicit GSkyRegionMap(const GSkyMap& map);
    explicit GSkyRegionMap(const GSkyRegion* region);
    GSkyRegionMap(const GSkyRegionMap& region);
    virtual ~GSkyRegionMap(void);

    // Operators
    GSkyRegionMap& operator=(const GSkyRegionMap& region);

    // Implemented pure virtual  methods
    void           clear(void);
    GSkyRegionMap* clone(void) const;
    std::string    classname(void) const;
    void           read(const std::string& line);
    std::string    write(void) const;
    bool           contains(const GSkyDir& dir) const;
    bool           contains(const GSkyRegion& reg) const;
    bool           overlaps(const GSkyRegion& reg) const;
    std::string    print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void                    load(const GFilename& filename);
    void                    map(const GSkyMap& map);
    const GSkyMap&          map(void) const;
    const std::vector<int>& nonzero_indices(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GSkyRegionMap& region);
    void free_members(void);
    void compute_solid_angle(void);
    void set_nonzero_indices(void);
    void set_region_circle(const GSkyRegionCircle* circle);
    void set_region_rectangle(const GSkyRegionRectangle* rect);

    // Protected members
    GSkyMap	          m_map;             //!< The region map
    std::vector<int>  m_nonzero_indices; //!< Vector of non-zero pixel indices
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GSkyRegionMap").
 ***************************************************************************/
inline
std::string GSkyRegionMap::classname(void) const
{
    return ("GSkyRegionMap");
}


/***********************************************************************//**
 * @brief Return sky map
 *
* @return region sky map.
 ***************************************************************************/
inline
const GSkyMap& GSkyRegionMap::map(void) const
{
    return (m_map);
}


/***********************************************************************//**
 * @brief Get non-zero index vector
 *
 * @return Reference to non-zero pixel indices vector.
 ***************************************************************************/
inline
const std::vector<int>& GSkyRegionMap::nonzero_indices(void) const
{
    return (m_nonzero_indices);
}


/***********************************************************************//**
 * @brief Return file name
 *
 * @return File name from which the region was loaded or into which
 *         the region was saved.
 *
 * Returns the file name from which the region was loaded or into
 * which the region was saved. The returned string will be empty if
 * no load() or save() method has been called before.
 ***************************************************************************/
/*
inline
const GFilename& GSkyRegionMap::filename(void) const
{
    return (m_filename);
}
*/

#endif /* GSKYREGIONMAP_HPP */
