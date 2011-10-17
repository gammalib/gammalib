/***************************************************************************
 *      GModelTemporalRegistry.hpp  -  Temporal model registry class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
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
 * @file GModelTemporalRegistry.hpp
 * @brief GModelTemporalRegistry class interface definition
 * @author J. Knodlseder
 */

#ifndef GMODELTEMPORALREGISTRY_HPP
#define GMODELTEMPORALREGISTRY_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GLog.hpp"
#include "GModelTemporal.hpp"


/***********************************************************************//**
 * @class GModelTemporalRegistry
 *
 * @brief Interface definition for the temporal model registry class.
 *
 * The registry class allows the registration of temporal models that are not
 * necessarily compiled into the GammaLib. It uses the static members
 * m_number, m_names, and m_models which are allocated globally to keep track
 * of temporal models that are available throughout all linked libraries. To
 * register a temporal model it is sufficient to add
 *  const GModelTemporalXXX      g_temporal_XXX_seed;
 *  const GModelTemporalRegistry g_temporal_XXX_registry(&g_spectral_XXX_seed);
 * at the top of the .cpp file of the temporal model. Here, XXX is a unique
 * name that describes the model.
 ***************************************************************************/
class GModelTemporalRegistry {

    // I/O friends
    friend std::ostream& operator<<(std::ostream& os, const GModelTemporalRegistry& registry);
    friend GLog&         operator<<(GLog& log, const GModelTemporalRegistry& registry);

public:
    // Constructors and destructors
    GModelTemporalRegistry(void);
    GModelTemporalRegistry(const GModelTemporal* model);
    GModelTemporalRegistry(const GModelTemporalRegistry& registry);
    virtual ~GModelTemporalRegistry(void);

    // Operators
    GModelTemporalRegistry& operator= (const GModelTemporalRegistry& registry);

    // Methods
    int             size(void) const { return m_number; }
    GModelTemporal* alloc(const std::string& type) const;
    std::string     name(const int& index) const;
    std::string     print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelTemporalRegistry& registry);
    void free_members(void);

private:
    // Pricate members
    static int                    m_number;   //!< Number of models in registry
    static std::string*           m_names;    //!< Model names
    static const GModelTemporal** m_models;   //!< Pointer to seed models
};

#endif /* GMODELTEMPORALREGISTRY_HPP */
