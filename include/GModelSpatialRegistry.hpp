/***************************************************************************
 *        GModelSpatialRegistry.hpp - Spatial model registry class         *
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
 * @file GModelSpatialRegistry.hpp
 * @brief Spatial model registry class definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPATIALREGISTRY_HPP
#define GMODELSPATIALREGISTRY_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GRegistry.hpp"
#include "GModelSpatial.hpp"


/***********************************************************************//**
 * @class GModelSpatialRegistry
 *
 * @brief Interface definition for the spatial model registry class
 *
 * The registry class allows the registration of spatial models that are not
 * necessarily compiled into the GammaLib. It uses the static members
 * m_number, m_names, and m_models which are allocated globally to keep track
 * of spatial models that are available throughout all linked libraries. To
 * register a spatial model it is sufficient to add
 *
 *     const GModelSpatialXXX      g_spatial_XXX_seed;
 *     const GModelSpatialRegistry g_spatial_XXX_registry(&g_spatial_XXX_seed);
 *
 * at the top of the .cpp file of the spatial model. Here, XXX is a unique
 * name that describes the model.
 ***************************************************************************/
class GModelSpatialRegistry : public GRegistry {

public:
    // Constructors and destructors
    GModelSpatialRegistry(void);
    GModelSpatialRegistry(const GModelSpatial* model);
    GModelSpatialRegistry(const GModelSpatialRegistry& registry);
    virtual ~GModelSpatialRegistry(void);

    // Operators
    GModelSpatialRegistry& operator=(const GModelSpatialRegistry& registry);

    // Methods
    int            size(void) const;
    GModelSpatial* alloc(const std::string& name) const;
    std::string    name(const int& index) const;
    std::string    print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpatialRegistry& registry);
    void free_members(void);

private:
    // Pricate members
    static int                   m_number;   //!< Number of models in registry
    static std::string*          m_names;    //!< Model names
    static const GModelSpatial** m_models;   //!< Pointer to seed models
};


/***********************************************************************//**
 * @brief Return number of registered models
 *
 * @return Number of registered models.
 *
 * Returns the number of registered model.
 ***************************************************************************/
inline
int GModelSpatialRegistry::size(void) const
{
    return m_number;
}

#endif /* GMODELSPATIALREGISTRY_HPP */
