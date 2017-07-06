/***************************************************************************
 *                  GSkyRegions.i - Sky region container class             *
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
 * @file GSkyRegions.i
 * @brief Sky regions container class interface file
 * @author Pierrick Martin
 */

%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSkyRegions.hpp"
%}

/***********************************************************************//**
 * @class GSkyRegions
 *
 * @brief Sky region container class
 ***************************************************************************/
class GSkyRegions : public GContainer {
public:
    // Constructors and destructors
    GSkyRegions(void);
    GSkyRegions(const GSkyRegions& regions);
	explicit GSkyRegions(const GFilename& filename);
    virtual ~GSkyRegions(void);

    // Methods
    void              clear(void);
    GSkyRegions*      clone(void) const;
    std::string       classname(void) const;
    GSkyRegion*       at(const int& index);
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
};

/***********************************************************************//**
 * @brief GSkyRegions class extension
 ***************************************************************************/
%extend GSkyRegions {
    GSkyRegion* __getitem__(const int& index) {
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", index, self->size());
        }
    }
    void __setitem__(const int& index, const GSkyRegion& val) {
        if (index>=0 && index < self->size()) {
            self->set(index, val);
            return;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", index, self->size());
        }
    }
    GSkyRegions copy() {
        return (*self);
    }
};
