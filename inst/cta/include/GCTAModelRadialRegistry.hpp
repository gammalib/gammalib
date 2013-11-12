/***************************************************************************
 *      GCTAModelRadialRegistry.hpp - CTA Radial model registry class      *
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
 * @file GCTAModelRadialRegistry.hpp
 * @brief CTA radial model registry class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAMODELRADIALREGISTRY_HPP
#define GCTAMODELRADIALREGISTRY_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GRegistry.hpp"
#include "GCTAModelRadial.hpp"


/***********************************************************************//**
 * @class GCTAModelRadialRegistry
 *
 * @brief Interface definition for the CTA radial model registry class
 *
 * The registry class allows the registration of radial models for CTA that
 * are not necessarily compiled into the GammaLib. It uses the static members
 * m_number, m_names, and m_models which are allocated globally to keep track
 * of radial models that are available throughout all linked libraries. To
 * register a radial model it is sufficient to add
 *
 *     const GCTAModelRadialXXX      g_cta_radial_XXX_seed;
 *     const GCTAModelRadialRegistry g_cta_radial_XXX_registry(&g_cta_radial_XXX_seed);
 *
 * at the top of the .cpp file of the radial model. Here, XXX is a unique
 * name that describes the model.
 ***************************************************************************/
class GCTAModelRadialRegistry : public GRegistry {

public:
    // Constructors and destructors
    GCTAModelRadialRegistry(void);
    GCTAModelRadialRegistry(const GCTAModelRadial* model);
    GCTAModelRadialRegistry(const GCTAModelRadialRegistry& registry);
    virtual ~GCTAModelRadialRegistry(void);

    // Operators
    GCTAModelRadialRegistry& operator= (const GCTAModelRadialRegistry& registry);

    // Methods
    int              size(void) const;
    GCTAModelRadial* alloc(const std::string& name) const;
    std::string      name(const int& index) const;
    std::string      print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAModelRadialRegistry& registry);
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
    static GRegistryPointer<const GCTAModelRadial*>& models() {
        static GRegistryPointer<const GCTAModelRadial*> m_models;
        return m_models;
    };
};


/***********************************************************************//**
 * @brief Return number of registered models
 *
 * @return Number of registered models.
 *
 * Returns the number of registered model.
 ***************************************************************************/
inline
int GCTAModelRadialRegistry::size(void) const
{
    return number();
}

#endif /* GCTAMODELRADIALREGISTRY_HPP */
