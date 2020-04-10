/***************************************************************************
 *             GResponseCacheKey.cpp - Response cache key class            *
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
 * @file GResponseCacheKey.cpp
 * @brief Response cache key class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GResponseCacheKey.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_SET                          "GResponseCacheKey::set(int&, float&)"
#define G_INSERT                    "GResponseCacheKey::insert(int&, float&)"
#define G_REMOVE                            "GResponseCacheKey::remove(int&)"

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
GResponseCacheKey::GResponseCacheKey(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Key name constructor
 *
 * @param[in] name Key name.
 ***************************************************************************/
GResponseCacheKey::GResponseCacheKey(const std::string& name)
{
    // Initialise class members
    init_members();

    // Set name
    this->name(name);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] key Response cache key.
 ***************************************************************************/
GResponseCacheKey::GResponseCacheKey(const GResponseCacheKey& key)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(key);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GResponseCacheKey::~GResponseCacheKey(void)
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
 * @param[in] key Response cache key.
 * @return Response cache key.
 ***************************************************************************/
GResponseCacheKey& GResponseCacheKey::operator=(const GResponseCacheKey& key)
{
    // Execute only if object is not identical
    if (this != &key) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(key);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Compare operator
 *
 * @param[in] key Response cache key.
 * @return True if instance is smaller than @p key.
 ***************************************************************************/
bool GResponseCacheKey::operator<(const GResponseCacheKey& key) const
{
    // Get key size
    int sz = size();

    // Compare key sizes
    bool smaller = sz < key.size();

    // If sizes are identical then continue witl comparison
    if (sz == key.size()) {

        // Compare key name
        smaller = m_name < key.m_name;

        // If names are identical then compare values. For small sizes we
        // explicitly unroll the loop to make sure that the code is fast.
        if (m_name == key.m_name) {
            if (sz == 1) {
                smaller = m_values[0] < key.m_values[0];
            }
            else if (sz == 2) {
                smaller = (m_values[0] < key.m_values[0]) ||
                          (m_values[1] < key.m_values[1]);
            }
            else if (sz == 3) {
                smaller = (m_values[0] < key.m_values[0]) ||
                          (m_values[1] < key.m_values[1]) ||
                          (m_values[2] < key.m_values[2]);
            }
            else if (sz == 4) {
                smaller = (m_values[0] < key.m_values[0]) ||
                          (m_values[1] < key.m_values[1]) ||
                          (m_values[2] < key.m_values[2]) ||
                          (m_values[3] < key.m_values[3]);
            }
            else if (sz == 5) {
                smaller = (m_values[0] < key.m_values[0]) ||
                          (m_values[1] < key.m_values[1]) ||
                          (m_values[2] < key.m_values[2]) ||
                          (m_values[3] < key.m_values[3]) ||
                          (m_values[4] < key.m_values[4]);
            }
            else {
                for (int i = 0; i < size(); ++i) {
                    smaller = m_values[i] < key.m_values[i];
                    if (smaller) {
                        break;
                    }
                }
            }
        } // endif: names were identical

    } // endif: sizes were identical

    // Return smaller flag
    return smaller;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear response cache
 ***************************************************************************/
void GResponseCacheKey::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 *
 * @return Pointer to deep copy of response cache key
 *
 * Makes a deep copy of the response cache key.
 ***************************************************************************/
GResponseCacheKey* GResponseCacheKey::clone(void) const
{
    return new GResponseCacheKey(*this);
}


/***********************************************************************//**
 * @brief Set key value
 *
 * @param[in] index Key value index [0,...,size()-1].
 * @param[in] value Key value.
 *
 * @exception GException::out_of_range
 *            Invalid key value index.
 *
 * Set key value in the container.
 ***************************************************************************/
void GResponseCacheKey::set(const int& index, const float& value)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_SET, "Key value index", index, size());
    }
    #endif

    // Assign key value
    m_values[index] = value;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Insert key value
 *
 * @param[in] index Key value index [0,...,size()-1].
 * @param[in] value Key value.
 *
 * @exception GException::out_of_range
 *            Invalid key value index.
 *
 * Inserts a key value before the value with the specified @p index.
 ***************************************************************************/
void GResponseCacheKey::insert(const int& index, const float& value)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_INSERT, "Key value index", index, size());
    }
    #endif

    // Inserts key value
    m_values.insert(m_values.begin()+index, value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Remove key value
 *
 * @param[in] index Key value index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Invalid key value index.
 *
 * Removes key value with the specified @p index.
 ***************************************************************************/
void GResponseCacheKey::remove(const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_REMOVE, "Key value index", index, size());
    }
    #endif

    // Erase key value
    m_values.erase(m_values.begin()+index);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print response cache key
 *
 * @param[in] chatter Chattiness.
 * @return String containing response cache key information.
 ***************************************************************************/
std::string GResponseCacheKey::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GResponseCacheKey ===");

        // Append total cache size
        result.append("\n"+gammalib::parformat("Key name"));
        result.append(m_name);
        result.append("\n"+gammalib::parformat("Number of key values"));
        result.append(gammalib::str(size()));
        for (int i = 0; i < size(); ++i) {
            result.append("\n"+gammalib::parformat("Key value "+gammalib::str(i+1)));
            result.append(gammalib::str(m_values[i]));
        }

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GResponseCacheKey::init_members(void)
{
    // Initialise members
    m_name.clear();
    m_values.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] key Response cache key.
 ***************************************************************************/
void GResponseCacheKey::copy_members(const GResponseCacheKey& key)
{
    // Copy members
    m_name   = key.m_name;
    m_values = key.m_values;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GResponseCacheKey::free_members(void)
{
    // Return
    return;
}
