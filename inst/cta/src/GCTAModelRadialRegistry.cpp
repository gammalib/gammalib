/***************************************************************************
 *      GCTAModelRadialRegistry.cpp - CTA Radial model registry class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2021 by Juergen Knoedlseder                         *
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
 * @file GCTAModelRadialRegistry.cpp
 * @brief CTA radial model registry class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GCTAModelRadialRegistry.hpp"
#include "GException.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_NAME                          "GCTAModelRadialRegistry::name(int&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
#define G_DEBUG_REGISTRY 0


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTAModelRadialRegistry::GCTAModelRadialRegistry(void)
{
    // Initialise members
    init_members();

    // Debug option: Show actual registry
    #if G_DEBUG_REGISTRY
    std::cout << "GCTAModelRadialRegistry(void): ";
    for (int i = 0; i < size(); ++i) {
        std::cout << "\"" << names()[i] << "\" ";
    }
    std::cout << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Model constructor
 *
 * @param[in] model Model.
 ***************************************************************************/
GCTAModelRadialRegistry::GCTAModelRadialRegistry(const GCTAModelRadial* model)
{
    // Initialise members
    init_members();

    // Debug option: Notify new registry
    #if G_DEBUG_REGISTRY
    std::cout << "GCTAModelRadialRegistry(const GCTAModelRadial*): ";
    std::cout << "add \"" << model->type() << "\" to registry." << std::endl;
    #endif

    // Allocate new registry
    std::string*            new_names  = new std::string[size()+1];
    const GCTAModelRadial** new_models = new const GCTAModelRadial*[size()+1];

    // Save old registry
    for (int i = 0; i < size(); ++i) {
        new_names[i]  = names()[i];
        new_models[i] = models()[i];
    }

    // Add new model to registry
    new_names[size()]  = model->type();
    new_models[size()] = model;

    // Set pointers on new registry
    names().assign(new_names);
    models().assign(new_models);

    // Increment number of models in registry
    number()++;

    // Debug option: Show actual registry
    #if G_DEBUG_REGISTRY
    std::cout << "GCTAModelRadialRegistry(const GCTAModelRadial*): ";
    for (int i = 0; i < size(); ++i) {
        std::cout << "\"" << names()[i] << "\" ";
    }
    std::cout << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] registry Registry.
 ***************************************************************************/
GCTAModelRadialRegistry::GCTAModelRadialRegistry(const GCTAModelRadialRegistry& registry)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(registry);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAModelRadialRegistry::~GCTAModelRadialRegistry(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] registry Registry.
 * @return Reference to registry.
 ***************************************************************************/
GCTAModelRadialRegistry& GCTAModelRadialRegistry::operator=(const GCTAModelRadialRegistry& registry)
{
    // Execute only if object is not identical
    if (this != &registry) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(registry);

    } // endif: object was not identical

    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Allocate radial model of given name
 *
 * @param[in] name Radial model name.
 *
 * Returns a pointer to a radial model instance of the specified name.
 * If the name has not been found in the registry, a NULL pointer is
 * returned.
 ***************************************************************************/
GCTAModelRadial* GCTAModelRadialRegistry::alloc(const std::string& name) const
{
    // Initialise radial model
    GCTAModelRadial* model = NULL;

    // Search for model in registry
    for (int i = 0; i < size(); ++i) {
        if (names()[i] == name) {
            model = models()[i]->clone();
            break;
        }
    }

    // Return radial model
    return model;
}


/***********************************************************************//**
 * @brief Returns model name
 *
 * @param[in] index Model index [0,...,size()-1].
 * @return Model name.
 *
 * @exception GException::out_of_range
 *            Model index is out of range.
 ***************************************************************************/
std::string GCTAModelRadialRegistry::name(const int& index) const
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_NAME, "Model index", index, size());
    }
    #endif

    // Return name
    return (names()[index]);
}


/***********************************************************************//**
 * @brief Print registry information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing registry information.
 ***************************************************************************/
std::string GCTAModelRadialRegistry::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAModelRadialRegistry ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of models"));
        result.append(gammalib::str(size()));

        // NORMAL: Append models
        if (chatter >= NORMAL) {
            for (int i = 0; i < size(); ++i) {
                result.append("\n"+gammalib::parformat(names()[i]));
                result.append(models()[i]->type());
            }
        }

    } // endif: chatter was not silent

    // Return result
    return result;
}

/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCTAModelRadialRegistry::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] registry Registry.
 ***************************************************************************/
void GCTAModelRadialRegistry::copy_members(const GCTAModelRadialRegistry& registry)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAModelRadialRegistry::free_members(void)
{
    // Return
    return;
}
