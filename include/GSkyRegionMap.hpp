/***************************************************************************
 *              GSkyRegionMap.hpp - sky map region class                *
 * ----------------------------------------------------------------------- *
 * copyright (C) 2013 by Michael Mayer                                     *
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
 * @brief Circular sky region class interface definition
 * @author Pierrick Martin and Anneli Schulz
 */

#ifndef GSKYREGIONMAP_HPP
#define GSKYREGIONMAP_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GMath.hpp"
#include "GSkyDir.hpp"
#include "GSkyRegion.hpp"


/***********************************************************************//**
 * @class GSkyRegionMap
 *
 * @brief Interface for a sky region in the form of a map
 *
 * This class provides an implementation for a sky region defined by a skymap.
 * The map is provided as a FITS file filled with 0 and 1 (or non-zero) values.
 *
 ***************************************************************************/
class GSkyRegionMap : public GSkyRegion {

public:
    // Constructors and destructors
    GSkyRegionMap(void);
    GSkyRegionMap(const GFits& file);
    GSkyRegionMap(const GFitsHDU& hdu);
	GSkyRegionMap(const GSkyRegionMap& map);
    virtual ~GSkyRegionMap(void);

    // Operators
    GSkyRegionMap& operator=(const GSkyRegionMap& map);

    // Implemented methods
    void              clear(void);
    GSkyRegionMap*    clone(void) const;
    void              load(const GFits& file);
    void              read(const GFitsHDU& hdu);   
    const GSkymap*    map(void) const;
	bool              contains(const GSkyDir& dir) const;
    bool              contains(const GSkyRegion& reg) const;
    bool              overlaps(const GSkyRegion& reg) const;
    std::string       print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GSkyRegionMap& region);
    void free_members(void);
    void compute_solid_angle(void);

    // Protected members
    GSkymap	  m_map;          //!< The map
    
};


/***********************************************************************//**
 * @brief Return the map
 *
 * @return GSkymap object.
 *
 * Returns the map.
 ***************************************************************************/
inline
double GSkyRegionMap::dec(void) const
{
    return (m_map);
}

#endif /* GSKYREGIONMAP_HPP */
