/***************************************************************************
 *             GCTAResponseCache.cpp - CTA response cache class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2018 by Juergen Knoedlseder                              *
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
 * @file GCTAResponseCache.cpp
 * @brief CTA response cache class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GCTAResponseCache.hpp"

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
GCTAResponseCache::GCTAResponseCache(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] cache CTA response cache.
 ***************************************************************************/
GCTAResponseCache::GCTAResponseCache(const GCTAResponseCache& cache)
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
GCTAResponseCache::~GCTAResponseCache(void)
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
 * @param[in] cache CTA response cache.
 * @return CTA response cache.
 ***************************************************************************/
GCTAResponseCache& GCTAResponseCache::operator=(const GCTAResponseCache& cache)
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
 * @brief Clear CTA response cache
 ***************************************************************************/
void GCTAResponseCache::clear(void)
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
int GCTAResponseCache::size(void) const
{
    // Initialize size
    int size = 0;

    // Compute size
    for (GCTAResponseCacheName::const_iterator it_name = m_cache.begin();
         it_name != m_cache.end(); ++it_name) {
        for (GCTAResponseCacheEreco::const_iterator it_ereco = it_name->second.begin();
             it_ereco != it_name->second.end(); ++it_ereco) {
            size += it_ereco->second.size();
        }
    }

    // Return size
    return size;
}


/***********************************************************************//**
 * @brief Return number of reconstructed energies in cache
 *
 * @return Number of reconstructed energies in cache
 *
 * Returns the number of reconstructed energies in the response cache.
 ***************************************************************************/
int GCTAResponseCache::nerecos(void) const
{
    // Initialize number of reconstructed energies
    int nerecos = 0;

    // Compute number of elements
    for (GCTAResponseCacheName::const_iterator it_name = m_cache.begin();
         it_name != m_cache.end(); ++it_name) {
        for (GCTAResponseCacheEreco::const_iterator it_ereco = it_name->second.begin();
             it_ereco != it_name->second.end(); ++it_ereco) {
            nerecos++;
        }
    }

    // Return number of reconstructed energies
    return nerecos;
}


/***********************************************************************//**
 * @brief Return number of true energies in cache
 *
 * @return Number of true energies in cache
 *
 * Returns the number of true energies in the response cache.
 ***************************************************************************/
int GCTAResponseCache::netrues(void) const
{
    // Initialize number of true energies
    int netrues = 0;

    // Compute number of elements
    for (GCTAResponseCacheName::const_iterator it_name = m_cache.begin();
         it_name != m_cache.end(); ++it_name) {
        for (GCTAResponseCacheEreco::const_iterator it_ereco = it_name->second.begin();
             it_ereco != it_name->second.end(); ++it_ereco) {
            netrues += it_ereco->second.size();
        }
    }

    // Return number of true energies
    return netrues;
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
void GCTAResponseCache::set(const std::string& name,
                            const GEnergy&     ereco,
                            const GEnergy&     etrue,
                            const double&      value)
{
    // Set cache value
    m_cache[name][ereco][etrue] = value;

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
 * Check if the cache contains a value for a specific set of parameters, and
 * optionally returns that value through the @p value pointer argument. If
 * @p value is set to NULL, no cached value is returned.
 ***************************************************************************/
bool GCTAResponseCache::contains(const std::string& name,
                                 const GEnergy&     ereco,
                                 const GEnergy&     etrue,
                                 double*            value) const
{
    // Initialise containment flag
    bool contains = false;

    // Search for name in cache
    GCTAResponseCacheName::const_iterator it_name = m_cache.find(name);
    if (it_name != m_cache.end()) {

        // Search for reconstructed energy in cache
        GCTAResponseCacheEreco::const_iterator it_ereco;
        it_ereco = it_name->second.find(ereco);
        if (it_ereco != it_name->second.end()) {

            // Search for true energy in cache
            GCTAResponseCacheEtrue::const_iterator it_etrue = it_ereco->second.find(etrue);
            if (it_etrue != it_ereco->second.end()) {
                contains = true;
                if (value != NULL) {
                    *value = it_etrue->second;
                }
            } // endif: true energy found

        } // endif: reconstructed energy found

    } // endif: name found

    // Return containment flag
    return contains;
}


/***********************************************************************//**
 * @brief Clone CTA response cache
 *
 * @return Pointer to deep copy of CTA response cache.
 ***************************************************************************/
GCTAResponseCache* GCTAResponseCache::clone(void) const
{
    return new GCTAResponseCache(*this);
}


/***********************************************************************//**
 * @brief Print CTA response cache
 *
 * @param[in] chatter Chattiness.
 * @return String containing CTA response cache information.
 *
 * @todo Implement method.
 ***************************************************************************/
std::string GCTAResponseCache::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAResponseCache ===");

        // Append total cache size
        result.append("\n"+gammalib::parformat("Number of cached values"));
        result.append(gammalib::str(m_cache.size()));

        // Append information
        for (GCTAResponseCacheName::const_iterator it_name = m_cache.begin();
             it_name != m_cache.end(); ++it_name) {

            // Append name
            result.append("\n"+gammalib::parformat("Name")+it_name->first);

            // Compute number of reconstructed and true energies
            int nereco = 0;
            int netrue = 0;
            for (GCTAResponseCacheEreco::const_iterator it_ereco = it_name->second.begin();
                 it_ereco != it_name->second.end(); ++it_ereco) {
                nereco++;
                netrue += it_ereco->second.size();
            }

            // Append number of reconstructed and true energies
            result.append("\n"+gammalib::parformat("Reconstructed energies"));
            result.append(gammalib::str(nereco));
            result.append("\n"+gammalib::parformat("True energies"));
            result.append(gammalib::str(netrue));

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
void GCTAResponseCache::init_members(void)
{
    // Initialise members
    m_cache.clear();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] cache CTA response cache.
 ***************************************************************************/
void GCTAResponseCache::copy_members(const GCTAResponseCache& cache)
{
    // Copy members
    m_cache = cache.m_cache;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAResponseCache::free_members(void)
{
    // Return
    return;
}
