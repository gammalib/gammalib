/***************************************************************************
 *        GModelSpatialRegistry.cpp - Spatial model registry class         *
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
 * @file GModelSpatialRegistry.cpp
 * @brief Spatial model registry class implemenation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GModelSpatialRegistry.hpp"
#include "GException.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_NAME                            "GModelSpatialRegistry::name(int&)"

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
GModelSpatialRegistry::GModelSpatialRegistry(void)
{
    // Initialise private members for clean destruction
    init_members();

    // Debug option: Show actual registry
    #if G_DEBUG_REGISTRY
    std::cout << "GModelSpatialRegistry(void): ";
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
 *
 * Construct registry by adding a model to the registry. This is the standard
 * constructor that is used to register a new model to GammaLib.
 ***************************************************************************/
GModelSpatialRegistry::GModelSpatialRegistry(const GModelSpatial* model)
{
    // Initialise private members for clean destruction
    init_members();

    // Debug option: Notify new registry
    #if G_DEBUG_REGISTRY
    std::cout << "GModelSpatialRegistry(const GModelSpatial*): ";
    std::cout << "add \"" << model->type() << "\" to registry." << std::endl;
    #endif

    // Allocate new registry
    std::string*          new_names  = new std::string[size()+1];
    const GModelSpatial** new_models = new const GModelSpatial*[size()+1];

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
    std::cout << "GModelSpatialRegistry(const GModelSpatial*): ";
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
GModelSpatialRegistry::GModelSpatialRegistry(const GModelSpatialRegistry& registry)
{
    // Initialise private members
    init_members();

    // Copy members
    copy_members(registry);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GModelSpatialRegistry::~GModelSpatialRegistry(void)
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
GModelSpatialRegistry& GModelSpatialRegistry::operator= (const GModelSpatialRegistry& registry)
{
    // Execute only if object is not identical
    if (this != &registry) {

        // Free members
        free_members();

        // Initialise private members
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
 * @brief Allocate spatial model of given name
 *
 * @param[in] name Model name.
 * @return Pointer to model (NULL if name is not registered).
 *
 * Returns a pointer to a void spatial model instance of the specified
 * name. If the name has not been found in the registry, a NULL pointer is
 * returned.
 ***************************************************************************/
GModelSpatial* GModelSpatialRegistry::alloc(const std::string& name) const
{
    // Initialise spatial model
    GModelSpatial* model = NULL;

    // Search for model in registry
    for (int i = 0; i < size(); ++i) {
        if (names()[i] == name) {
            model = models()[i]->clone();
            break;
        }
    }

    // Return spatial model
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
std::string GModelSpatialRegistry::name(const int& index) const
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_NAME, index, 0, size()-1);
    }
    #endif

    // Return name
    return (names()[index]);
}


/***********************************************************************//**
 * @brief Print registry information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return Registry content.
 ***************************************************************************/
std::string GModelSpatialRegistry::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpatialRegistry ===");

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
void GModelSpatialRegistry::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] registry Registry.
 ***************************************************************************/
void GModelSpatialRegistry::copy_members(const GModelSpatialRegistry& registry)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatialRegistry::free_members(void)
{
    // Return
    return;
}
