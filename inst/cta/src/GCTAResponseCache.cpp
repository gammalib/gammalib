/***************************************************************************
 *             GCTAResponseCache.cpp - CTA response cache class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2018-2019 by Juergen Knoedlseder                         *
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
#include "GCTAInstDir.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_NAME                                "GCTAResponseCache::name(int&)"

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
        for (GCTAResponseCacheRA::const_iterator it_ra =
             it_name->second.begin();
             it_ra != it_name->second.end(); ++it_ra) {
            for (GCTAResponseCacheDEC::const_iterator it_dec =
                 it_ra->second.begin();
                 it_dec != it_ra->second.end(); ++it_dec) {
                for (GCTAResponseCacheEreco::const_iterator it_ereco =
                     it_dec->second.begin();
                     it_ereco != it_dec->second.end(); ++it_ereco) {
                    size += it_ereco->second.size();
                }
            }
        }
    }

    // Return size
    return size;
}


/***********************************************************************//**
 * @brief Return number of sky directions in cache
 *
 * @return Number of sky directions in cache
 *
 * Returns the number of sky directions in the response cache.
 ***************************************************************************/
int GCTAResponseCache::ndirs(void) const
{
    // Initialize number of reconstructed energies
    int ndirs = 0;

    // Compute number of elements
    for (GCTAResponseCacheName::const_iterator it_name = m_cache.begin();
         it_name != m_cache.end(); ++it_name) {
        for (GCTAResponseCacheRA::const_iterator it_ra = it_name->second.begin();
             it_ra != it_name->second.end(); ++it_ra) {
            for (GCTAResponseCacheDEC::const_iterator it_dec = it_ra->second.begin();
                 it_dec != it_ra->second.end(); ++it_dec) {
                ndirs += it_dec->second.size();
            }
        }
    }

    // Return number of sky directions
    return ndirs;
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
        for (GCTAResponseCacheRA::const_iterator it_ra = it_name->second.begin();
             it_ra != it_name->second.end(); ++it_ra) {
            for (GCTAResponseCacheDEC::const_iterator it_dec = it_ra->second.begin();
                 it_dec != it_ra->second.end(); ++it_dec) {
                for (GCTAResponseCacheEreco::const_iterator it_ereco =
                     it_dec->second.begin();
                     it_ereco != it_dec->second.end(); ++it_ereco) {
                    nerecos++;
                }
            }
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
        for (GCTAResponseCacheRA::const_iterator it_ra = it_name->second.begin();
             it_ra != it_name->second.end(); ++it_ra) {
            for (GCTAResponseCacheDEC::const_iterator it_dec = it_ra->second.begin();
                 it_dec != it_ra->second.end(); ++it_dec) {
                for (GCTAResponseCacheEreco::const_iterator it_ereco =
                     it_dec->second.begin();
                     it_ereco != it_dec->second.end(); ++it_ereco) {
                    netrues += it_ereco->second.size();
                }
            }
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
 *
 * The method will set the Right Ascension and Declination of the cached
 * value to zero.
 ***************************************************************************/
void GCTAResponseCache::set(const std::string& name,
                            const GEnergy&     ereco,
                            const GEnergy&     etrue,
                            const double&      value)
{
    // Set cache value
    m_cache[name][0.0][0.0][ereco][etrue] = value;

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
void GCTAResponseCache::set(const std::string& name,
                            const GCTAInstDir& dir,
                            const GEnergy&     ereco,
                            const GEnergy&     etrue,
                            const double&      value)
{
    // Get Right Ascension and Declination of instrument direction
    const double& ra  = dir.dir().ra_deg();
    const double& dec = dir.dir().dec_deg();

    // Set cache value
    m_cache[name][ra][dec][ereco][etrue] = value;

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
 * The method assumes that the Right Ascension and Declination of the cached
 * value is zero.
 *
 * If the @p value pointer argument is not NULL, the method will return the
 * cached value through this argument in case that the value exists.
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

        // Search for Right Ascension in cache
        GCTAResponseCacheRA::const_iterator it_ra = it_name->second.find(0.0);
        if (it_ra != it_name->second.end()) {

            // Search for Declination in cache
            GCTAResponseCacheDEC::const_iterator it_dec = it_ra->second.find(0.0);
            if (it_dec != it_ra->second.end()) {

                // Search for reconstructed energy in cache
                GCTAResponseCacheEreco::const_iterator it_ereco =
                                        it_dec->second.find(ereco);
                if (it_ereco != it_dec->second.end()) {

                    // Search for true energy in cache
                    GCTAResponseCacheEtrue::const_iterator it_etrue =
                                            it_ereco->second.find(etrue);
                    if (it_etrue != it_ereco->second.end()) {
                        contains = true;
                        if (value != NULL) {
                            *value = it_etrue->second;
                        }
                    } // endif: true energy found

                } // endif: reconstructed energy found

            } // endif: Declination found

        } // endif: Right Ascension found

    } // endif: name found

    // Return containment flag
    return contains;
}


/***********************************************************************//**
 * @brief Check if cache contains a value for specific parameters
 *
 * @param[in] name Cache name.
 * @param[in] dir Sky direction.
 * @param[in] ereco Reconstructed energy.
 * @param[in] etrue True energy.
 * @param[out] value Pointer to cached value (only if found).
 * @return True if cached value was found, otherwise false.
 *
 * Check if the cache contains a value for a given @p name, instrument
 * direction @p dir,reconstructed energy @p ereco, and true energy @p etrue.
 *
 * If the @p value pointer argument is not NULL, the method will return the
 * cached value through this argument in case that the value exists.
 ***************************************************************************/
bool GCTAResponseCache::contains(const std::string& name,
                                 const GCTAInstDir& dir,
                                 const GEnergy&     ereco,
                                 const GEnergy&     etrue,
                                 double*            value) const
{
    // Initialise containment flag
    bool contains = false;

    // Get Right Ascension and Declination of instrument direction
    const double& ra  = dir.dir().ra_deg();
    const double& dec = dir.dir().dec_deg();

    // Search for name in cache
    GCTAResponseCacheName::const_iterator it_name = m_cache.find(name);
    if (it_name != m_cache.end()) {

        // Search for Right Ascension in cache
        GCTAResponseCacheRA::const_iterator it_ra = it_name->second.find(ra);
        if (it_ra != it_name->second.end()) {

            // Search for Declination in cache
            GCTAResponseCacheDEC::const_iterator it_dec = it_ra->second.find(dec);
            if (it_dec != it_ra->second.end()) {

                // Search for reconstructed energy in cache
                GCTAResponseCacheEreco::const_iterator it_ereco =
                                        it_dec->second.find(ereco);
                if (it_ereco != it_dec->second.end()) {

                    // Search for true energy in cache
                    GCTAResponseCacheEtrue::const_iterator it_etrue =
                                            it_ereco->second.find(etrue);
                    if (it_etrue != it_ereco->second.end()) {
                        contains = true;
                        if (value != NULL) {
                            *value = it_etrue->second;
                        }
                    } // endif: true energy found

                } // endif: reconstructed energy found

            } // endif: Declination found

        } // endif: Right Ascension found

    } // endif: name found

    // Return containment flag
    return contains;
}


/***********************************************************************//**
 * @brief Return name of cache element
 *
 * @param[in] index Cache element number (0,...,nnames()-1).
 * @return Name of cache element.
 *
 * Returns name of cache element with specified @p index.
 ***************************************************************************/
std::string GCTAResponseCache::name(const int& index) const
{
    // Raise an exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_NAME, index, size());
    }

    // Initialise name
    std::string name;

    // Initialise counter
    int counter = 0;

    // Loop over names and extract the name that corresponds to the index
    for (GCTAResponseCacheName::const_iterator it_name = m_cache.begin();
         it_name != m_cache.end(); ++it_name, ++counter) {
        if (counter == index) {
            name = it_name->first;
            break;
        }
    }

    // Return name
    return name;
}


/***********************************************************************//**
 * @brief Remove cache values for a given name
 *
 * @param[in] name Cache name.
 *
 * Removes the cache values for a given @p name. If the @p name does not
 * exist in the cache the method does nothing.
 ***************************************************************************/
void GCTAResponseCache::remove(const std::string& name)
{
    // Search for name in cache
    GCTAResponseCacheName::const_iterator it_name = m_cache.find(name);
    if (it_name != m_cache.end()) {
        m_cache.erase(it_name);
    }

    // Return
    return;
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
        result.append(gammalib::str(size()));

        // Append information
        for (GCTAResponseCacheName::const_iterator it_name = m_cache.begin();
             it_name != m_cache.end(); ++it_name) {

            // Append name
            result.append("\n"+gammalib::parformat("Name")+it_name->first);

            // Compute number of sky directions, reconstructed and true energies
            int ndir   = 0;
            int nereco = 0;
            int netrue = 0;
            for (GCTAResponseCacheRA::const_iterator it_ra =
                 it_name->second.begin();
                 it_ra != it_name->second.end(); ++it_ra) {
                for (GCTAResponseCacheDEC::const_iterator it_dec =
                     it_ra->second.begin();
                     it_dec != it_ra->second.end(); ++it_dec) {
                    ndir++;
                    for (GCTAResponseCacheEreco::const_iterator it_ereco =
                         it_dec->second.begin();
                         it_ereco != it_dec->second.end(); ++it_ereco) {
                        nereco++;
                        netrue += it_ereco->second.size();
                    }
                }
            }

            // Append number of reconstructed and true energies
            result.append("\n"+gammalib::parformat("Instrument directions"));
            result.append(gammalib::str(ndir));
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
