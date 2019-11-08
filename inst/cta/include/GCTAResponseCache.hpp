/***************************************************************************
 *             GCTAResponseCache.hpp - CTA response cache class            *
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
 * @file GCTAResponseCache.hpp
 * @brief CTA response cache class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTARESPONSECACHE_HPP
#define GCTARESPONSECACHE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <map>
#include "GBase.hpp"
#include "GEnergy.hpp"

/* __ Forward declarations _______________________________________________ */
class GCTAInstDir;

/* __ Constants __________________________________________________________ */

/* __ Typedefs ___________________________________________________________ */
typedef std::map<double,double>                        GCTAResponseCacheElement;
typedef std::map<std::string,GCTAResponseCacheElement> GCTAResponseCacheName;


/***********************************************************************//**
 * @class GCTAResponseCache
 *
 * @brief CTA response cache class
 *
 * The class implements a cache for the Instrument Response Function values
 * so that the values do not need to be recomputed each time but can be
 * fetched from the cache.
 ***************************************************************************/
class GCTAResponseCache : public GBase {

public:
    // Constructors and destructors
    GCTAResponseCache(void);
    GCTAResponseCache(const GCTAResponseCache& cache);
    virtual ~GCTAResponseCache(void);

    // Operators
    GCTAResponseCache& operator=(const GCTAResponseCache& cache);

    // Methods
    void               clear(void);
    GCTAResponseCache* clone(void) const;
    std::string        classname(void) const;
    int                size(void) const;
    bool               is_empty(void) const;
    void               set(const std::string& name,
                           const GEnergy&     ereco,
                           const GEnergy&     etrue,
                           const double&      value);
    void               set(const std::string& name,
                           const GCTAInstDir& dir,
                           const GEnergy&     ereco,
                           const GEnergy&     etrue,
                           const double&      value);
    void               remove(const std::string& name);
    bool               contains(const std::string& name,
                                const GEnergy&     ereco,
                                const GEnergy&     etrue,
                                double*            value = NULL) const;
    bool               contains(const std::string& name,
                                const GCTAInstDir& dir,
                                const GEnergy&     ereco,
                                const GEnergy&     etrue,
                                double*            value = NULL) const;
    std::string        print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void   init_members(void);
    void   copy_members(const GCTAResponseCache& cache);
    void   free_members(void);
    double encode(const GEnergy&     ereco,
                  const GEnergy&     etrue) const;
    double encode(const GCTAInstDir& dir,
                  const GEnergy&     ereco,
                  const GEnergy&     etrue) const;

    // Protected members
    GCTAResponseCacheName m_cache;
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTAResponseCache").
 ***************************************************************************/
inline
std::string GCTAResponseCache::classname(void) const
{
    return ("GCTAResponseCache");
}


/***********************************************************************//**
 * @brief Checks whether the cache is empty
 *
 * @return True if cache is empty, false otherwise.
 *
 * Checks whether the response cache is empty.
 ***************************************************************************/
inline
bool GCTAResponseCache::is_empty(void) const
{
    return (m_cache.empty());
}


/***********************************************************************//**
 * @brief Remove cache for source
 *
 * @param[in] name Source name.
 *
 * Remove cache for source with @p name.
 ***************************************************************************/
inline
void GCTAResponseCache::remove(const std::string& name)
{
    m_cache.erase(name);
    return;
}

#endif /* GCTARESPONSECACHE_HPP */
