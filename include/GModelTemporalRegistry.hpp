/***************************************************************************
 *       GModelTemporalRegistry.hpp - Temporal model registry class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2016 by Juergen Knoedlseder                         *
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
 * @brief Temporal model registry class definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELTEMPORALREGISTRY_HPP
#define GMODELTEMPORALREGISTRY_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GRegistry.hpp"

/* __ Forward declarations _______________________________________________ */
class GXmlElement;
class GModelTemporal;


/***********************************************************************//**
 * @class GModelTemporalRegistry
 *
 * @brief Interface definition for the temporal model registry class.
 *
 * The registry class allows the registration of temporal models that are not
 * necessarily compiled into the GammaLib. It uses the static members
 * m_number and m_models which are allocated globally to keep track of
 * temporal models that are available throughout all linked libraries. To
 * register a temporal model it is sufficient to add
 *
 *     const GModelTemporalXXX      g_temporal_XXX_seed;
 *     const GModelTemporalRegistry g_temporal_XXX_registry(&g_spectral_XXX_seed);
 *
 * at the top of the .cpp file of the temporal model. Here, XXX is a unique
 * name that describes the model.
 ***************************************************************************/
class GModelTemporalRegistry : public GRegistry {

public:
    // Constructors and destructors
    GModelTemporalRegistry(void);
    GModelTemporalRegistry(const GModelTemporal* model);
    GModelTemporalRegistry(const GModelTemporalRegistry& registry);
    virtual ~GModelTemporalRegistry(void);

    // Operators
    GModelTemporalRegistry& operator=(const GModelTemporalRegistry& registry);

    // Methods
    std::string     classname(void) const;
    int             size(void) const;
    GModelTemporal* alloc(const GXmlElement& xml) const;
    std::string     name(const int& index) const;
    std::string     print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelTemporalRegistry& registry);
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
    // Pointer to seed models
    static GRegistryPointer<const GModelTemporal*>& models() {
        static GRegistryPointer<const GModelTemporal*> m_models;
        return m_models;
    };
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelTemporalRegistry").
 ***************************************************************************/
inline
std::string GModelTemporalRegistry::classname(void) const
{
    return ("GModelTemporalRegistry");
}


/***********************************************************************//**
 * @brief Return number of registered models
 *
 * @return Number of registered models.
 *
 * Returns the number of registered model.
 ***************************************************************************/
inline
int GModelTemporalRegistry::size(void) const
{
    return number();
}

#endif /* GMODELTEMPORALREGISTRY_HPP */
