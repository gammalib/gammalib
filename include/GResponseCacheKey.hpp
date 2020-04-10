/***************************************************************************
 *            GResponseCacheKey.hpp - Response cache key class             *
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
 * @file GResponseCacheKey.hpp
 * @brief Response cache key class definition
 * @author Juergen Knoedlseder
 */

#ifndef GRESPONSECACHEKEY_HPP
#define GRESPONSECACHEKEY_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GContainer.hpp"

/* __ Forward declarations _______________________________________________ */

/* __ Constants __________________________________________________________ */

/* __ Typedefs ___________________________________________________________ */


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
    void               insert(const int& index, const float& value);
    void               remove(const int& index);
    void               reserve(const int& num);
    void               name(const std::string& name);
    const std::string& name(void) const;
    std::string        print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void   init_members(void);
    void   copy_members(const GResponseCacheKey& key);
    void   free_members(void);

    // Protected members
    std::string        m_name;    //!< Key name (e.g. model name, observation)
    std::vector<float> m_values;  //!< Key values
};


/***********************************************************************//**
 * @brief Return reference to key value
 *
 * @param[in] index Key index.
 *
 * Returns a reference to a key value.
 ***************************************************************************/
float& GResponseCacheKey::operator[](const int& index)
{
    return m_values[index];
}


/***********************************************************************//**
 * @brief Return reference to key value (const version)
 *
 * @param[in] index Key index.
 *
 * Returns a reference to a key value.
 ***************************************************************************/
const float& GResponseCacheKey::operator[](const int& index) const
{
    return m_values[index];
}


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GResponseCacheKey").
 ***************************************************************************/
inline
std::string GResponseCacheKey::classname(void) const
{
    return ("GResponseCacheKey");
}


/***********************************************************************//**
 * @brief Checks whether the key has no values
 *
 * @return True if key has no values, false otherwise.
 *
 * Checks whether the key has no values.
 ***************************************************************************/
inline
bool GResponseCacheKey::is_empty(void) const
{
    return (m_values.empty());
}


/***********************************************************************//**
 * @brief Return number of key values
 *
 * @return Number of key values.
 *
 * Returns the number of key values.
 ***************************************************************************/
inline
int GResponseCacheKey::size(void) const
{
    return (m_values.size());
}


/***********************************************************************//**
 * @brief Append key value
 *
 * @param[in] value Key value.
 *
 * Append @p value to key values.
 ***************************************************************************/
inline
void GResponseCacheKey::append(const float& value)
{
    m_values.push_back(value);
    return;
}


/***********************************************************************//**
 * @brief Reserves space for key values
 *
 * @param[in] num Number of key values
 *
 * Reserves space for @p num key values.
 ***************************************************************************/
inline
void GResponseCacheKey::reserve(const int& num)
{
    m_values.reserve(num);
    return;
}


/***********************************************************************//**
 * @brief Set key name
 *
 * @param[in] name Key name
 *
 * Sets the key name.
 ***************************************************************************/
inline
void GResponseCacheKey::name(const std::string& name)
{
    m_name = name;
    return;
}


/***********************************************************************//**
 * @brief Return key name
 *
 * @returns Key name
 *
 * Returns the key name.
 ***************************************************************************/
inline
const std::string& GResponseCacheKey::name(void) const
{
    return m_name;
}

#endif /* GRESPONSECACHEKEY_HPP */
