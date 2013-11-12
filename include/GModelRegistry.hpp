/***************************************************************************
 *                GModelRegistry.hpp - Model registry class                *
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
 * @file GModelRegistry.hpp
 * @brief Model registry class definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELREGISTRY_HPP
#define GMODELREGISTRY_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GRegistry.hpp"
#include "GModel.hpp"


/***********************************************************************//**
 * @class GModelRegistry
 *
 * @brief Interface definition for the model registry class
 *
 * The registry class allows the registration of models that are not
 * necessarily compiled into the GammaLib. It uses the static members
 * m_number, m_names, and m_models which are allocated globally to keep track
 * of models that are available throughout all linked libraries. To register
 * a model it is sufficient to add
 *
 *     const GModelXXX      g_XXX_seed;
 *     const GModelRegistry g_XXX_registry(&g_XXX_seed);
 *
 * at the top of the .cpp file of the model. Here, XXX is a unique name that
 * describes the model.
 ***************************************************************************/
class GModelRegistry : public GRegistry {

public:
    // Constructors and destructors
    GModelRegistry(void);
    GModelRegistry(const GModel* model);
    GModelRegistry(const GModelRegistry& registry);
    virtual ~GModelRegistry(void);

    // Operators
    GModelRegistry& operator=(const GModelRegistry& registry);

    // Methods
    int         size(void) const;
    GModel*     alloc(const std::string& name) const;
    std::string name(const int& index) const;
    std::string print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelRegistry& registry);
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
    static GRegistryPointer<const GModel*>& models() {
        static GRegistryPointer<const GModel*> m_models;
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
int GModelRegistry::size(void) const
{
    return number();
}

#endif /* GMODELREGISTRY_HPP */
