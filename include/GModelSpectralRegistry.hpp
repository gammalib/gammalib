/***************************************************************************
 *       GModelSpectralRegistry.hpp - Spectral model registry class        *
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
 * @file GModelSpectralRegistry.hpp
 * @brief Spectral model registry class definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPECTRALREGISTRY_HPP
#define GMODELSPECTRALREGISTRY_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GRegistry.hpp"

/* __ Forward declarations _______________________________________________ */
class GXmlElement;
class GModelSpectral;


/***********************************************************************//**
 * @class GModelSpectralRegistry
 *
 * @brief Interface definition for the spectral model registry class
 *
 * The registry class allows the registration of spectral models that are not
 * necessarily compiled into the GammaLib. It uses the static members
 * m_number and m_models which are allocated globally to keep track of
 * spectral models that are available throughout all linked libraries. To
 * register a spectral model it is sufficient to add
 *
 *     const GModelSpectralXXX      g_spectral_XXX_seed;
 *     const GModelSpectralRegistry g_spectral_XXX_registry(&g_spectral_XXX_seed);
 *
 * at the top of the .cpp file of the spectral model. Here, XXX is a unique
 * name that describes the model.
 ***************************************************************************/
class GModelSpectralRegistry : public GRegistry {

public:
    // Constructors and destructors
    GModelSpectralRegistry(void);
    GModelSpectralRegistry(const GModelSpectral* model);
    GModelSpectralRegistry(const GModelSpectralRegistry& registry);
    virtual ~GModelSpectralRegistry(void);

    // Operators
    GModelSpectralRegistry& operator=(const GModelSpectralRegistry& registry);

    // Methods
    std::string     classname(void) const;
    int             size(void) const;
    GModelSpectral* alloc(const GXmlElement& xml) const;
    std::string     name(const int& index) const;
    std::string     print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralRegistry& registry);
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
    static GRegistryPointer<const GModelSpectral*>& models() {
        static GRegistryPointer<const GModelSpectral*> m_models;
        return m_models;
    };
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpectralRegistry").
 ***************************************************************************/
inline
std::string GModelSpectralRegistry::classname(void) const
{
    return ("GModelSpectralRegistry");
}


/***********************************************************************//**
 * @brief Return number of registered models
 *
 * @return Number of registered models.
 *
 * Returns the number of registered model.
 ***************************************************************************/
inline
int GModelSpectralRegistry::size(void) const
{
    return number();
}

#endif /* GMODELSPECTRALREGISTRY_HPP */
