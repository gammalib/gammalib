/***************************************************************************
 *          GResponseVectorCache.hpp - Response vector cache class         *
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
 * @file GResponseVectorCache.hpp
 * @brief Response vector cache class definition
 * @author Juergen Knoedlseder
 */

#ifndef GRESPONSEVECTORCACHE_HPP
#define GRESPONSEVECTORCACHE_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GBase.hpp"

/* __ Forward declarations _______________________________________________ */
class GVector;
class GFilename;
class GFitsTable;

/* __ Constants __________________________________________________________ */

/* __ Typedefs ___________________________________________________________ */


/***********************************************************************//**
 * @class GResponseVectorCache
 *
 * @brief Response vector cache class
 *
 * The class implements a vector cache for the Instrument Response Function
 * values so that the values do not need to be recomputed each time but can
 * be fetched from the cache.
 ***************************************************************************/
class GResponseVectorCache : public GBase {

public:
    // Constructors and destructors
    GResponseVectorCache(void);
    GResponseVectorCache(const GResponseVectorCache& cache);
    virtual ~GResponseVectorCache(void);

    // Operators
    GResponseVectorCache& operator=(const GResponseVectorCache& cache);

    // Methods
    void                  clear(void);
    GResponseVectorCache* clone(void) const;
    std::string           classname(void) const;
    bool                  is_empty(void) const;
    int                   size(void) const;
    void                  set(const std::string& cache_id,
                              const GVector&     vector);
    void                  remove(const std::string& cache_id);
    bool                  contains(const std::string& cache_id,
                                   GVector*           irfs = NULL) const;
    void                  load(const GFilename& filename);
    void                  save(const GFilename& filename,
                               const bool& clobber = false) const;
    void                  read(const GFitsTable& table);
    std::string           print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GResponseVectorCache& cache);
    void free_members(void);
    int  find_cache(const std::string& cache_id) const;

    // Protected members
    std::vector<std::string> m_cache_ids;
    std::vector<int>         m_cache_entries;
    std::vector<double*>     m_cache_values;
    std::vector<int*>        m_cache_indices;
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GResponseVectorCache").
 ***************************************************************************/
inline
std::string GResponseVectorCache::classname(void) const
{
    return ("GResponseVectorCache");
}


/***********************************************************************//**
 * @brief Checks whether the cache is empty
 *
 * @return True if cache is empty, false otherwise.
 *
 * Checks whether the response cache is empty.
 ***************************************************************************/
inline
bool GResponseVectorCache::is_empty(void) const
{
    return (m_cache_ids.empty());
}


/***********************************************************************//**
 * @brief Returns size of vector chache
 *
 * @return Size of vector cache.
 *
 * Returns the number of vectors that are stored in the vector cache.
 ***************************************************************************/
inline
int GResponseVectorCache::size(void) const
{
    return ((int)m_cache_ids.size());
}

#endif /* GRESPONSEVECTORCACHE_HPP */
