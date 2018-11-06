/***************************************************************************
 *             GCTAResponseCache.hpp - CTA response cache class            *
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
typedef std::map<GEnergy,double>                  GCTAResponseCacheEtrue;
typedef std::map<GEnergy,GCTAResponseCacheEtrue>  GCTAResponseCacheEreco;
typedef std::map<double,GCTAResponseCacheEreco>   GCTAResponseCacheDEC;
typedef std::map<double,GCTAResponseCacheDEC>     GCTAResponseCacheRA;
typedef std::map<std::string,GCTAResponseCacheRA> GCTAResponseCacheName;


/***********************************************************************//**
 * @class GCTAResponseCache
 *
 * @brief CTA response cache class
 *
 * @todo Add class description.
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
    int                ndirs(void) const;
    int                nerecos(void) const;
    int                netrues(void) const;
    void               set(const std::string& name,
                           const GEnergy&     ereco,
                           const GEnergy&     etrue,
                           const double&      value);
    void               set(const std::string& name,
                           const GCTAInstDir& dir,
                           const GEnergy&     ereco,
                           const GEnergy&     etrue,
                           const double&      value);
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
    void init_members(void);
    void copy_members(const GCTAResponseCache& cache);
    void free_members(void);

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

#endif /* GCTARESPONSECACHE_HPP */
