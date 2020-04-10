/***************************************************************************
 *             GResponseCacheKey.i - Response cache key class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GResponseCacheKey.i
 * @brief Response cache key class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GResponseCacheKey.hpp"
%}


/***********************************************************************//**
 * @class GResponseCacheKey
 *
 * @brief Response cache key class
 *
 * The class implements the key for the Instrument Response Function cache.
 * The key aggregates all relevant attributes to describe a cached value,
 * such as model name, true and measured event energies and instrument
 * direction.
 ***************************************************************************/
class GResponseCacheKey : public GContainer {

public:
    // Constructors and destructors
    GResponseCacheKey(void);
    explicit GResponseCacheKey(const std::string& name);
    GResponseCacheKey(const GResponseCacheKey& key);
    virtual ~GResponseCacheKey(void);

    // Operators
    GResponseCacheKey& operator=(const GResponseCacheKey& key);
    float&             operator[](const int& index);
    const float&       operator[](const int& index) const;
    bool               operator<(const GResponseCacheKey& key) const;

    // Methods
    void               clear(void);
    GResponseCacheKey* clone(void) const;
    std::string        classname(void) const;
    bool               is_empty(void) const;
    int                size(void) const;
    void               set(const int& index, const float& value);
    void               append(const float& value);
    void               insert(const int& index, float& value);
    void               remove(const int& index);
    void               reserve(const int& num);
    void               name(const std::string& name);
    const std::string& name(void) const;
};


/***********************************************************************//**
 * @brief GResponseCacheKey class extension
 ***************************************************************************/
%extend GResponseCacheKey {
    GResponseCacheKey copy() {
        return (*self);
    }
};
