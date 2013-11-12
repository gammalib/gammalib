/***************************************************************************
 *          GObservationRegistry.hpp - Observation registry class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
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
 * @file GObservationRegistry.hpp
 * @brief Observation registry class definition
 * @author Juergen Knoedlseder
 */

#ifndef GOBSERVATIONREGISTRY_HPP
#define GOBSERVATIONREGISTRY_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GRegistry.hpp"
#include "GObservation.hpp"


/***********************************************************************//**
 * @class GObservationRegistry
 *
 * @brief Interface definition for the observation registry class
 *
 * The registry class allows the registration of observations that are not
 * necessarily compiled into the GammaLib. It uses the static members
 * m_number, m_names, and m_obs which are allocated globally to keep track
 * of observations that are available throughout all linked libraries. To
 * register an observation it is sufficient to add
 *
 *     const GXXXObservation      g_obs_XXX_seed;
 *     const GObservationRegistry g_obs_XXX_registry(&g_obs_XXX_seed);
 *
 * at the top of the .cpp file of the observation. Here, XXX is a unique
 * name that describes the instrument for which the observation class is
 * implemented.
 ***************************************************************************/
class GObservationRegistry : public GRegistry {

public:
    // Constructors and destructors
    GObservationRegistry(void);
    GObservationRegistry(const GObservation* obs);
    GObservationRegistry(const GObservationRegistry& registry);
    virtual ~GObservationRegistry(void);

    // Operators
    GObservationRegistry& operator= (const GObservationRegistry& registry);

    // Methods
    int           size(void) const;
    GObservation* alloc(const std::string& name) const;
    std::string   name(const int& index) const;
    std::string   print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GObservationRegistry& registry);
    void free_members(void);

private:
    // Private members (the private members have been implement as static
    // methods to avoid the static initialization order fiasco of static
    // members; using static methods we follow the "construct on first use
    // idiom")
    // Number of models in registry
    static int& number() {
        static int m_number = 0;
        return m_number;
    };
    // Model names
    static GRegistryPointer<std::string>& names() {
        static GRegistryPointer<std::string> m_names;
        return m_names;
    };
    // Pointer to seed models
    static GRegistryPointer<const GObservation*>& obs() {
        static GRegistryPointer<const GObservation*> m_obs;
        return m_obs;
    };
};


/***********************************************************************//**
 * @brief Return number of registered observations
 *
 * @return Number of registered observations.
 *
 * Returns the number of registered observations.
 ***************************************************************************/
inline
int GObservationRegistry::size(void) const
{
    return number();
}

#endif /* GOBSERVATIONREGISTRY_HPP */
