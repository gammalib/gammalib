/***************************************************************************
 *                GResponseCache.cpp - Response cache class                *
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
 * @file GResponseCache.cpp
 * @brief Response cache class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GResponseCache.hpp"
#include "GInstDir.hpp"

/* __ Method name definitions ____________________________________________ */

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
GResponseCache::GResponseCache(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] cache Response cache.
 ***************************************************************************/
GResponseCache::GResponseCache(const GResponseCache& cache)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(cache);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GResponseCache::~GResponseCache(void)
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
 * @param[in] cache Response cache.
 * @return Response cache.
 ***************************************************************************/
GResponseCache& GResponseCache::operator=(const GResponseCache& cache)
{
    // Execute only if object is not identical
    if (this != &cache) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(cache);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear response cache
 ***************************************************************************/
void GResponseCache::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return number of elements in cache
 *
 * @return Number of elements in cache
 *
 * Returns the number of elements in the response cache.
 ***************************************************************************/
int GResponseCache::size(void) const
{
    // Initialize size
    int size = 0;

    // Compute size
    for (GResponseCacheName::const_iterator it_name = m_cache.begin();
         it_name != m_cache.end(); ++it_name) {
        for (GResponseCacheEnergy::const_iterator it_energy =
             it_name->second.begin();
             it_energy != it_name->second.end(); ++it_energy) {
                size += it_energy->second.size();
/*
            for (GResponseCacheDir::const_iterator it_dir =
                 it_energy->second.begin();
                 it_dir != it_energy->second.end(); ++it_dir) {
                size += it_dir->second.size();
            }
*/
        }
    }

    // Return size
    return size;
}


/***********************************************************************//**
 * @brief Return number of energies in cache
 *
 * @return Number of energies in cache
 *
 * Returns the number of energies in the response cache.
 ***************************************************************************/
int GResponseCache::nenergies(void) const
{
    // Initialize number of energies
    int nenergies = 0;

    // Compute number of elements
    for (GResponseCacheName::const_iterator it_name = m_cache.begin();
         it_name != m_cache.end(); ++it_name) {
        for (GResponseCacheEnergy::const_iterator it_energy =
             it_name->second.begin();
             it_energy != it_name->second.end(); ++it_energy) {
            nenergies += it_energy->second.size();
        }
    }

    // Return number of energies
    return nenergies;
}


/***********************************************************************//**
 * @brief Set cache value
 *
 * @param[in] name Cache name.
 * @param[in] ereco Reconstructed energy.
 * @param[in] etrue True energy.
 * @param[in] value Cache value.
 *
 * Set cache value for a given @p name, reconstructed energy @p ereco, and
 * true energy @p etrue.
 ***************************************************************************/
void GResponseCache::set(const std::string& name,
                         const GEnergy&     ereco,
                         const GEnergy&     etrue,
                         const double&      value)
{
    // Get energy identifier
    u_int64_t energy = hash(ereco, etrue);

    // Set cache value
    m_cache[name][energy][0] = value;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set cache value
 *
 * @param[in] name Cache name.
 * @param[in] dir Instrument direction.
 * @param[in] ereco Reconstructed energy.
 * @param[in] etrue True energy.
 * @param[in] value Cache value.
 *
 * Set cache value for a given @p name, instrument direction @p dir,
 * reconstructed energy @p ereco, and true energy @p etrue.
 ***************************************************************************/
void GResponseCache::set(const std::string& name,
                         const GInstDir&    dir,
                         const GEnergy&     ereco,
                         const GEnergy&     etrue,
                         const double&      value)
{
    // Get energy and direction identifier
    u_int64_t energy  = hash(ereco, etrue);
    u_int64_t instdir = dir.hash();

    // Set cache value
    m_cache[name][energy][instdir] = value;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Check if cache contains a value for specific parameters
 *
 * @param[in] name Cache name.
 * @param[in] ereco Reconstructed energy.
 * @param[in] etrue True energy.
 * @param[out] value Pointer to cached value (only if found).
 * @return True if cached value was found, otherwise false.
 *
 * Check if the cache contains a value for a given @p name, reconstructed
 * energy @p ereco, and true energy @p etrue.
 *
 * If the @p value pointer argument is not NULL, the method will return the
 * cached value through this argument in case that the value exists.
 ***************************************************************************/
bool GResponseCache::contains(const std::string& name,
                              const GEnergy&     ereco,
                              const GEnergy&     etrue,
                              double*            value) const
{
    // Initialise containment flag
    bool contains = false;

    // Get energy identifier
    u_int64_t energy  = hash(ereco, etrue);
    u_int64_t instdir = 0;

    // Search for name in cache
    GResponseCacheName::const_iterator it_name = m_cache.find(name);
    if (it_name != m_cache.end()) {

        // Search for energy in cache
        GResponseCacheEnergy::const_iterator it_energy =
                                  it_name->second.find(energy);
        if (it_energy != it_name->second.end()) {

            // Search for direction in cache
            GResponseCacheDir::const_iterator it_dir =
                                  it_energy->second.find(instdir);
            if (it_dir != it_energy->second.end()) {
                contains = true;
                if (value != NULL) {
                    *value = it_dir->second;
                }
            } // endif: direction found

        } // endif: energy found

    } // endif: name found

    // Return containment flag
    return contains;
}


/***********************************************************************//**
 * @brief Check if cache contains a value for specific parameters
 *
 * @param[in] name Cache name.
 * @param[in] dir Instrument direction.
 * @param[in] ereco Reconstructed energy.
 * @param[in] etrue True energy.
 * @param[out] value Pointer to cached value (only if found).
 * @return True if cached value was found, otherwise false.
 *
 * Check if the cache contains a value for a given @p name, instrument
 * direction @p dir, reconstructed energy @p ereco, and true energy @p etrue.
 *
 * If the @p value pointer argument is not NULL, the method will return the
 * cached value through this argument in case that the value exists.
 ***************************************************************************/
bool GResponseCache::contains(const std::string& name,
                              const GInstDir&    dir,
                              const GEnergy&     ereco,
                              const GEnergy&     etrue,
                              double*            value) const
{
    // Initialise containment flag
    bool contains = false;

    // Get energy identifier
    u_int64_t energy  = hash(ereco, etrue);
    u_int64_t instdir = dir.hash();

    // Search for name in cache
    GResponseCacheName::const_iterator it_name = m_cache.find(name);
    if (it_name != m_cache.end()) {

        // Search for energy in cache
        GResponseCacheEnergy::const_iterator it_energy =
                                  it_name->second.find(energy);
        if (it_energy != it_name->second.end()) {

            // Search for direction in cache
            GResponseCacheDir::const_iterator it_dir =
                                  it_energy->second.find(instdir);
            if (it_dir != it_energy->second.end()) {
                contains = true;
                if (value != NULL) {
                    *value = it_dir->second;
                }
            } // endif: direction found

        } // endif: energy found

    } // endif: name found

    // Return containment flag
    return contains;
}


/***********************************************************************//**
 * @brief Clone response cache
 *
 * @return Pointer to deep copy of response cache.
 ***************************************************************************/
GResponseCache* GResponseCache::clone(void) const
{
    return new GResponseCache(*this);
}


/***********************************************************************//**
 * @brief Print response cache
 *
 * @param[in] chatter Chattiness.
 * @return String containing response cache information.
 ***************************************************************************/
std::string GResponseCache::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GResponseCache ===");

        // Append total cache size
        result.append("\n"+gammalib::parformat("Number of cached values"));
        result.append(gammalib::str(size()));

        // Append information
        for (GResponseCacheName::const_iterator it_name = m_cache.begin();
             it_name != m_cache.end(); ++it_name) {

            // Append name
            result.append("\n"+gammalib::parformat("Name")+it_name->first);

            // Compute number of sky directions and energies
            int ndirs     = 0;
            int nenergies = 0;
            for (GResponseCacheEnergy::const_iterator it_energy =
                 it_name->second.begin();
                 it_energy != it_name->second.end(); ++it_energy) {
                nenergies++;
                for (GResponseCacheDir::const_iterator it_dir =
                     it_energy->second.begin();
                     it_dir != it_energy->second.end(); ++it_dir) {
                    ndirs++;
                }
            }

            // Append number of reconstructed and true energies
            result.append("\n"+gammalib::parformat("Number of energies"));
            result.append(gammalib::str(nenergies));
            result.append("\n"+gammalib::parformat("Number of directions"));
            result.append(gammalib::str(ndirs));

        } // endfor: looped over names

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
void GResponseCache::init_members(void)
{
    // Initialise members
    m_cache.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] cache Response cache.
 ***************************************************************************/
void GResponseCache::copy_members(const GResponseCache& cache)
{
    // Copy members
    m_cache = cache.m_cache;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GResponseCache::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Encode reconstructued and true energy into hash value
 *
 * @param[in] ereco Reconstructed energy.
 * @param[in] etrue True energy.
 * @return Encoded reconstructued and true energy
 *
 * Encodes the reconstructued and true energy in a 64 Bit unsigned integer
 * value. The energy are converted into single precision floating point
 * values and mapped on a 64 Bit unsigned integer value. This provides unique
 * hash values for all instruments up to the precision of floating point
 * singles.
 ***************************************************************************/
u_int64_t GResponseCache::hash(const GEnergy&  ereco,
                               const GEnergy&  etrue) const
{
    // Allocate static array to store the two energies as floats
    static float buffer[2];

    // Store the two energies as floats
    buffer[0] = float(ereco.MeV());
    buffer[1] = float(etrue.MeV());

    // Map the floats to an unsigned 64 Bit integer
    u_int64_t hash; std::memcpy(&hash, &buffer, sizeof hash);

    // Return encoded value
    return hash;
}
