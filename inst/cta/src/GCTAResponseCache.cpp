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
        size += it_name->second.size();
    }

    // Return size
    return size;
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
    // Get element identifier
    double element = encode(ereco, etrue);

    // Set cache value
    m_cache[name][element] = value;

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
    // Get element identifier
    double element = encode(dir, ereco, etrue);

    // Set cache value
    m_cache[name][element] = value;

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
bool GCTAResponseCache::contains(const std::string& name,
                                 const GEnergy&     ereco,
                                 const GEnergy&     etrue,
                                 double*            value) const
{
    // Initialise containment flag
    bool contains = false;

    // Get element identifier
    double element = encode(ereco, etrue);

    // Search for name in cache
    GCTAResponseCacheName::const_iterator it_name = m_cache.find(name);
    if (it_name != m_cache.end()) {

        // Search for element in cache
        GCTAResponseCacheElement::const_iterator it_element =
                                  it_name->second.find(element);
        if (it_element != it_name->second.end()) {
            contains = true;
            if (value != NULL) {
                *value = it_element->second;
            }
        } // endif: element found

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

    // Get element identifier
    double element = encode(dir, ereco, etrue);

    // Search for name in cache
    GCTAResponseCacheName::const_iterator it_name = m_cache.find(name);
    if (it_name != m_cache.end()) {

        // Search for element in cache
        GCTAResponseCacheElement::const_iterator it_element =
                                  it_name->second.find(element);
        if (it_element != it_name->second.end()) {
            contains = true;
            if (value != NULL) {
                *value = it_element->second;
            }
        } // endif: element found

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
        result.append(gammalib::str(size()));

        // Append information
        for (GCTAResponseCacheName::const_iterator it_name = m_cache.begin();
             it_name != m_cache.end(); ++it_name) {

            // Append name
            result.append("\n"+gammalib::parformat("Name")+it_name->first);

            // Compute number of elements
            int nelement = it_name->second.size();

            // Append number of elements
            result.append("\n"+gammalib::parformat("Elements"));
            result.append(gammalib::str(nelement));

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


/***********************************************************************//**
 * @brief Encode reconstructued and true energy
 *
 * @param[in] ereco Reconstructed energy.
 * @param[in] etrue True energy.
 * @return Encoded reconstructued and true energy
 *
 * Encodes the reconstructued and true energy in a single double precision
 * value. The encoding is done using the following formula:
 *
 * \f[
 *    {\\ encode} = E_{\rm reco} \times 10^2 + E_{\rm true}
 * \f]
 *
 * where
 * - \f$E_{\rm reco}\f$ is the reconstructued energy in TeV, and
 * - \f$E_{\rm true}\f$ is the true energy in TeV.
 ***************************************************************************/
double GCTAResponseCache::encode(const GEnergy& ereco,
                                 const GEnergy& etrue) const
{
    // Construct encoded instrument direction and reconstructued energy
    double encoded = ereco.TeV() * 1.0e2 + etrue.TeV();

    // Return encoded value
    return encoded;
}


/***********************************************************************//**
 * @brief Encode instrument direction, reconstructued and true energy
 *
 * @param[in] dir Instrument direction.
 * @param[in] ereco Reconstructed energy.
 * @param[in] etrue True energy.
 * @return Encoded instrument direction, reconstructued and true energy
 *
 * Encodes the instrument direction, reconstructued and true energy in a
 * single double precision value. The encoding is done using the following
 * formula:
 *
 * \f[
 *    {\\ encode} = \alpha \times 10^8 + \delta \times 10^5 +
 *                  E_{\rm reco} \times 10^2 + E_{\rm true}
 * \f]
 *
 * where
 * - \f$\alpha\f$ is the Right Ascension of the instrument direction in
 *   degrees,
 * - \f$\delta\f$ is the Declination of the instrument direction in degrees,
 * - \f$E_{\rm reco}\f$ is the reconstructued energy in TeV, and
 * - \f$E_{\rm true}\f$ is the true energy in TeV.
 ***************************************************************************/
double GCTAResponseCache::encode(const GCTAInstDir& dir,
                                 const GEnergy&     ereco,
                                 const GEnergy&     etrue) const
{
    // Construct encoded instrument direction and reconstructued energy
    double encoded = dir.dir().ra_deg()  * 1.0e8 +
                     dir.dir().dec_deg() * 1.0e5 +
                     ereco.TeV()         * 1.0e2 +
                     etrue.TeV();

    // Return encoded value
    return encoded;
}
