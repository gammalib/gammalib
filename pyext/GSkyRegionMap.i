/***************************************************************************
 *                  GSkyRegionMap.i - Sky region map class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Pierrick Martin                                  *
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
 * @file GSkyRegionMap.i
 * @brief Sky region map class interface definition
 * @author Pierrick Martin
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSkyRegionMap.hpp"
%}


/***********************************************************************//**
 * @class GSkyRegionMap
 *
 * @brief Interface for the sky region map class
 ***************************************************************************/
class GSkyRegionMap : public GSkyRegion {
public:
    // Constructors and destructors
    GSkyRegionMap(void);
    GSkyRegionMap(const GFilename& filename);
    GSkyRegionMap(const GSkyMap& map); 
    GSkyRegionMap(const GSkyRegion* region);
	GSkyRegionMap(const GSkyRegionMap& region);
    virtual ~GSkyRegionMap(void);

    // Implemented methods
    void              clear(void);
    GSkyRegionMap*    clone(void) const;
    std::string       classname(void) const;
    void              read(const std::string& line);
    std::string       write(void) const;
    bool              contains(const GSkyDir& dir) const;
    bool              contains(const GSkyRegion& reg) const;
    bool              overlaps(const GSkyRegion& reg) const;
    
    // Other methods
    void                    load(const GFilename& filename);
    void                    map(const GSkyMap& map);
    const GSkyMap&          map(void) const;
    const GFilename&        filename(void) const;
    const std::vector<int>& nonzero_indices(void) const;
};


/***********************************************************************//**
 * @brief GSkyRegionMap class extension
 ***************************************************************************/
%extend GSkyRegionMap {
    GSkyRegionMap copy() {
        return (*self);
    }
};
