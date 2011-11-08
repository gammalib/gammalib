/***************************************************************************
 *         GObservationRegistry.hpp  -  Observation registry class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Juergen Knoedlseder                              *
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
 * @brief GObservationRegistry class interface definition
 * @author J. Knoedlseder
 */

#ifndef GOBSERVATIONREGISTRY_HPP
#define GOBSERVATIONREGISTRY_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GLog.hpp"
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
 *  const GXXXObservation      g_obs_XXX_seed;
 *  const GObservationRegistry g_obs_XXX_registry(&g_obs_XXX_seed);
 * at the top of the .cpp file of the observation. Here, XXX is a unique
 * name that describes the instrument for which the observation class is
 * implemented.
 ***************************************************************************/
class GObservationRegistry {

    // I/O friends
    friend std::ostream& operator<<(std::ostream& os, const GObservationRegistry& registry);
    friend GLog&         operator<<(GLog& log,        const GObservationRegistry& registry);

public:
    // Constructors and destructors
    GObservationRegistry(void);
    GObservationRegistry(const GObservation* obs);
    GObservationRegistry(const GObservationRegistry& registry);
    virtual ~GObservationRegistry(void);

    // Operators
    GObservationRegistry& operator= (const GObservationRegistry& registry);

    // Methods
    int           size(void) const { return m_number; }
    GObservation* alloc(const std::string& instrument) const;
    std::string   instrument(const int& index) const;
    std::string   print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GObservationRegistry& registry);
    void free_members(void);

private:
    // Private members
    static int                  m_number;  //!< Number of observations in registry
    static std::string*         m_names;   //!< Instrument names
    static const GObservation** m_obs;     //!< Pointer to seed observations
};

#endif /* GOBSERVATIONREGISTRY_HPP */
