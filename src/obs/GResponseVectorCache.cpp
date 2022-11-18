/***************************************************************************
 *         GResponseVectorCache.cpp - Response vector cache class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020-2022 by Juergen Knoedlseder                         *
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
 * @file GResponseVectorCache.cpp
 * @brief Response vector cache class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GVector.hpp"
#include "GFilename.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableStringCol.hpp"
#include "GFitsTableLongCol.hpp"
#include "GFitsTableDoubleCol.hpp"
#include "GResponseVectorCache.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_VECTOR_CACHE                        //!< Debug vector cache


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GResponseVectorCache::GResponseVectorCache(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] cache Response vector cache.
 ***************************************************************************/
GResponseVectorCache::GResponseVectorCache(const GResponseVectorCache& cache)
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
GResponseVectorCache::~GResponseVectorCache(void)
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
 * @param[in] cache Response vector cache.
 * @return Response vector cache.
 ***************************************************************************/
GResponseVectorCache& GResponseVectorCache::operator=(const GResponseVectorCache& cache)
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
 * @brief Clear response vector cache
 ***************************************************************************/
void GResponseVectorCache::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone response cache
 *
 * @return Pointer to deep copy of response cache.
 ***************************************************************************/
GResponseVectorCache* GResponseVectorCache::clone(void) const
{
    return new GResponseVectorCache(*this);
}


/***********************************************************************//**
 * @brief Set cache value
 *
 * @param[in] cache_id Cache identifier.
 * @param[in] vector Cache vector.
 *
 * Set cache vector for a given @p cache_id.
 ***************************************************************************/
void GResponseVectorCache::set(const std::string& cache_id,
                               const GVector&     vector)
{
    // Debug vector cache
    #if defined(G_DEBUG_VECTOR_CACHE)
    std::cout << "GResponseVectorCache::set(";
    std::cout << cache_id << "," << vector.size() << "):";
    #endif

    // Initialise pointers to values and indices
    double* values  = NULL;
    int*    indices = NULL;

    // Find cache index
    int index = find_cache(cache_id);

    // Determine size of vector cache
    int entries = vector.non_zeros();

    // Debug vector cache
    #if defined(G_DEBUG_VECTOR_CACHE)
    if (index != -1) {
        std::cout << " entry found at index " << index;
        std::cout << " with " << m_cache_entries[index] << " elements.";
        std::cout << " Require " << entries << " elements.";
        std::cout << std::endl;
    }
    else {
        std::cout << " no entry found.";
        std::cout << " Require " << entries << " elements.";
        std::cout << std::endl;
    }
    #endif

    // If cache identifier does not yet exist then create new cache entry
    // and allocate memory to hold cache values and indices
    if (index == -1) {
        if (entries > 0) {
            values  = new double[entries];
            indices = new int[entries];
        }
        m_cache_ids.push_back(cache_id);
        m_cache_entries.push_back(entries);
        m_cache_values.push_back(values);
        m_cache_indices.push_back(indices);
    }

    // ... otherwise update cache existing entry by allocating new memory
    // to hold cache values and indices in case that the number of entries
    // have changed
    else {
        if (m_cache_entries[index] != entries) {
            if (entries > 0) {
                values  = new double[entries];
                indices = new int[entries];
            }
            if (m_cache_values[index]  != NULL) delete [] m_cache_values[index];
            if (m_cache_indices[index] != NULL) delete [] m_cache_indices[index];
            m_cache_entries[index] = entries;
            m_cache_values[index]  = values;
            m_cache_indices[index] = indices;
        }
        else {
            values  = m_cache_values[index];
            indices = m_cache_indices[index];
        }
    }

    // Set cache vector
    for (int i = 0; i < vector.size(); ++i) {
        if (vector[i] != 0.0) {
            *values++  = vector[i];
            *indices++ = i;
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Remove cache
 *
 * @param[in] cache_id Cache identifier.
 *
 * Remove cache for a given @p cache_id.
 ***************************************************************************/
void GResponseVectorCache::remove(const std::string& cache_id)
{
    // Find cache index
    int index = find_cache(cache_id);

    // If cache index is valid then remove corresponding cache
    if (index != -1) {

        // Free memory
        if (m_cache_values[index]  != NULL) delete [] m_cache_values[index];
        if (m_cache_indices[index] != NULL) delete [] m_cache_indices[index];

        // Erase entry
        m_cache_ids.erase(m_cache_ids.begin() + index);
        m_cache_entries.erase(m_cache_entries.begin() + index);
        m_cache_values.erase(m_cache_values.begin() + index);
        m_cache_indices.erase(m_cache_indices.begin() + index);

    } // endif: cache index was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Check if cache contains a value for specific parameters
 *
 * @param[in] cache_id Cache identifier.
 * @param[in,out] vector Pointer to cached vector (only if found).
 * @return True if cached value was found, otherwise false.
 *
 * Check if the cache contains a value for a given @p name, reconstructed
 * energy @p ereco, and true energy @p etrue.
 *
 * If the @p value pointer argument is not NULL, the method will return the
 * cached value through this argument in case that the value exists.
 ***************************************************************************/
bool GResponseVectorCache::contains(const std::string& cache_id,
                                    GVector*           vector) const
{
    // Debug vector cache
    #if defined(G_DEBUG_VECTOR_CACHE)
    std::cout << "GResponseVectorCache::contains(" << cache_id;
    if (vector == NULL) {
        std::cout << ", NULL):";
    }
    else {
        std::cout << "," << vector->size() << "):";
    }
    #endif

    // Find cache index
    int index = find_cache(cache_id);

    // Set contains flag
    bool contains = (index != -1);

    // Debug vector cache
    #if defined(G_DEBUG_VECTOR_CACHE)
    if (contains) {
        std::cout << " entry found at index " << index;
        std::cout << " with " << m_cache_entries[index] << " elements.";
        std::cout << std::endl;
    }
    else {
        std::cout << " no entry found." << std::endl;
    }
    #endif

    // If cache contains the identifier and a pointer to a vector was
    // provided then fill the vector
    if (contains && vector != NULL) {

        // Get pointer to values and indices
        double* values  = m_cache_values[index];
        int*    indices = m_cache_indices[index];

        // Fill vector
        for (int i = 0; i < m_cache_entries[index]; ++i, ++values, ++indices) {
            (*vector)[*indices] = *values;
        }

    } // endif: filled vector

    // Return containment flag
    return contains;
}


/***********************************************************************//**
 * @brief Load response vector cache from FITS file
 *
 * @param[in] filename FITS file name.
 *
 * Loads the response vector cache from a FITS file. All binary tables in
 * the FITS file are assumed to be response vector cache entries.
 ***************************************************************************/
void GResponseVectorCache::load(const GFilename& filename)
{
    // Free members
    free_members();

    // Initialise attributes
    init_members();

    // Open FITS file
    GFits fits(filename);

    // Loop over all extensions
    for (int i = 0; i < fits.size(); ++i) {

        // Handle only binary table extensions
        if (fits[i]->exttype() == GFitsHDU::HT_BIN_TABLE) {
            read(*fits.table(i));
        }

    } // endfor: looped over all extensions

    // Close FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save the response vector cache into FITS file
 *
 * @param[in] filename FITS file name.
 * @param[in] clobber Overwrite an response vector cache file?
 *
 * Saves the response vector cache into a FITS file. If a file with the given
 * @p filename does not yet exist it will be created, otherwise the method
 * opens the existing file. The response vector cache can only be appended
 * to an existing file if the @p clobber flag is set to "true" (otherwise
 * an exception is thrown).
 *
 * The method will append all cache entries as binary FITS tables to the
 * FITS file. The table extension names are defined by the cache
 * identifiers. All cache entries with identifiers exceeding 80 characters
 * will be skipped to avoid truncation of cache identifiers.
 ***************************************************************************/
void GResponseVectorCache::save(const GFilename&   filename,
                                const bool&        clobber) const
{
    // Open or create FITS file (without extension name)
    GFits fits(filename.url(), true);

    // Loop over the cache entries
    for (int i = 0; i < size(); ++i) {

        // Set extension name
        std::string extname = m_cache_ids[i];

        // Skip entries with too long identifiers
        if (extname.length() > 80) {
            continue;
        }

        // Get number of cache entries
        int entries = m_cache_entries[i];

        // Create response vector cache columns
        GFitsTableDoubleCol col_values("VALUES",   entries);
        GFitsTableLongCol   col_indices("INDICES", entries);

        // Write cache values and indices in columns
        double* values  = m_cache_values[i];
        int*    indices = m_cache_indices[i];
        for (int k = 0; k < entries; ++k, ++values, ++indices) {
            col_values(k)  = *values;
            col_indices(k) = *indices;
        }

        // Create binary table
        GFitsBinTable table(entries);
        table.append(col_values);
        table.append(col_indices);
        table.extname(extname);

        // Write mandatory keywords
        table.card("HDUCLASS", "OGIP",         "Format conforms to OGIP standard");
        table.card("HDUCLAS1", "RESPONSE",     "Extension contains response data");
        table.card("HDUCLAS2", "VECTOR_CACHE", "Extension contains a response vector cache");
        table.card("HDUVERS",  "1.0.0",        "Version of the file format");

        // If the FITS object contains already an extension with the same
        // name then remove now this extension
        if (fits.contains(extname)) {
            fits.remove(extname);
        }

        // Append response vector cache to FITS file
        fits.append(table);


    } // endfor: looped over cache entries

    // Save to file
    fits.save(clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read response vector cache from FITS table
 *
 * @param[in] table FITS table.
 *
 * Reads response vector cache from a FITS table. There is one cache entry
 * in a FITS table, with the table extension name specifying the cache
 * identifier. Only tables containing the "HDUCLAS2" keyword set to
 * "VECTOR_CACHE" will be considered by the method, all other tables will
 * be simply skipped.
 ***************************************************************************/
void GResponseVectorCache::read(const GFitsTable& table)
{
    // Continue only if table contains a response vector cache
    if (table.has_card("HDUCLAS2") && table.string("HDUCLAS2") == "VECTOR_CACHE") {

        // Get number of rows in FITS table
        int entries = table.nrows();

        // Continue only if there are rows in the FITS table
        if (entries > 0) {

            // Set cache ID from extension name
            std::string cache_id = table.extname();

            // Get column pointers
            const GFitsTableCol* col_values  = table["VALUES"];
            const GFitsTableCol* col_indices = table["INDICES"];

            // Allocate memory for values and indices
            double* values  = new double[entries];
            int*    indices = new int[entries];

            // Put information and pointers into cache
            m_cache_ids.push_back(cache_id);
            m_cache_entries.push_back(entries);
            m_cache_values.push_back(values);
            m_cache_indices.push_back(indices);

            // Extact values and indices
            for (int i = 0; i < entries; ++i) {
                *values++  = col_values->real(i);
                *indices++ = col_indices->integer(i);
            }

        } // endif: there were entries in table

    } // endif: table contained a response vector cache

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print response cache
 *
 * @param[in] chatter Chattiness.
 * @return String containing response cache information.
 ***************************************************************************/
std::string GResponseVectorCache::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GResponseVectorCache ===");

        // Append total cache size
        result.append("\n"+gammalib::parformat("Number of cached vectors"));
        result.append(gammalib::str(size()));

        // Append cache information
        for (int i = 0; i < size(); ++i) {
            result.append("\n"+gammalib::parformat("Identifier"));
            result.append(m_cache_ids[i]);
            result.append(" ("+gammalib::str(m_cache_entries[i])+" values");
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
void GResponseVectorCache::init_members(void)
{
    // Initialise members
    m_cache_ids.clear();
    m_cache_entries.clear();
    m_cache_values.clear();
    m_cache_indices.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] cache Response cache.
 ***************************************************************************/
void GResponseVectorCache::copy_members(const GResponseVectorCache& cache)
{
    // Copy members
    m_cache_ids     = cache.m_cache_ids;
    m_cache_entries = cache.m_cache_entries;

    // Copy values and indices
    for (int i = 0; i < size(); ++i) {
        double* values  = NULL;
        int*    indices = NULL;
        if (m_cache_entries[i] > 0) {
            values  = new double[m_cache_entries[i]];
            indices = new int[m_cache_entries[i]];
            for (int k = 0; k < m_cache_entries[i]; ++k) {
                values[k]  = cache.m_cache_values[i][k];
                indices[k] = cache.m_cache_indices[i][k];
            }
        }
        m_cache_values.push_back(values);
        m_cache_indices.push_back(indices);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GResponseVectorCache::free_members(void)
{
    // Loop over entries
    for (int i = 0; i < size(); ++i) {
        if (m_cache_values[i]  != NULL) delete [] m_cache_values[i];
        if (m_cache_indices[i] != NULL) delete [] m_cache_indices[i];
        m_cache_values[i]  = NULL;
        m_cache_indices[i] = NULL;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Find cache
 *
 * @param[in] cache_id Cache identifier.
 * @return Cache index (-1 if no cache was found)
 *
 * Find cache for a given @p cache_id. If no cache was found the method
 * returns -1.
 ***************************************************************************/
int GResponseVectorCache::find_cache(const std::string& cache_id) const
{
    // Initialise cache index
    int index = -1;

    // Loop over cache
    for (int i = 0; i < m_cache_ids.size(); ++i) {
        if (m_cache_ids[i] == cache_id) {
            index    = i;
            break;
        }
    }

    // Return index
    return index;
}


