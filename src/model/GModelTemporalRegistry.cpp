/***************************************************************************
 *       GModelTemporalRegistry.cpp - Temporal model registry class        *
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
 * @file GModelTemporalRegistry.cpp
 * @brief Temporal model registry class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GModelTemporalRegistry.hpp"
#include "GException.hpp"
#include "GTools.hpp"

/* __ Static members _____________________________________________________ */
int                                     GModelTemporalRegistry::m_number(0);
//std::string*           GModelTemporalRegistry::m_names(0);
//const GModelTemporal** GModelTemporalRegistry::m_models(0);
GRegistryPointer<std::string>           GModelTemporalRegistry::m_names(0);
GRegistryPointer<const GModelTemporal*> GModelTemporalRegistry::m_models(0);

/* __ Method name definitions ____________________________________________ */
#define G_NAME                           "GModelTemporalRegistry::name(int&)"

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
GModelTemporalRegistry::GModelTemporalRegistry(void)
{
    // Initialise private members for clean destruction
    init_members();

    // Debug option: Show actual registry
    #if G_DEBUG_REGISTRY
    std::cout << "GModelTemporalRegistry(void): ";
    for (int i = 0; i < m_number; ++i) {
        std::cout << "\"" << m_names[i] << "\" ";
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
GModelTemporalRegistry::GModelTemporalRegistry(const GModelTemporal* model)
{
    // Initialise private members for clean destruction
    init_members();

    // Debug option: Notify new registry
    #if G_DEBUG_REGISTRY
    std::cout << "GModelTemporalRegistry(const GModelTemporal*): ";
    std::cout << "add \"" << model->type() << "\" to registry." << std::endl;
    #endif
    
    // Allocate new old registry
    std::string*           new_names  = new std::string[m_number+1];
    const GModelTemporal** new_models = new const GModelTemporal*[m_number+1];

    // Save old registry
    for (int i = 0; i < m_number; ++i) {
        new_names[i]  = m_names[i];
        new_models[i] = m_models[i];
    }

    // Add new model to registry
    new_names[m_number]  = model->type();
    new_models[m_number] = model;

    // Delete old registry
    //if (m_names  != NULL) delete [] m_names;
    //if (m_models != NULL) delete [] m_models;

    // Set pointers on new registry
    //m_names  = new_names;
    //m_models = new_models;
    m_names.assign(new_names);
    m_models.assign(new_models);

    // Increment number of models in registry
    m_number++;

    // Debug option: Show actual registry
    #if G_DEBUG_REGISTRY
    std::cout << "GModelTemporalRegistry(const GModelTemporal*): ";
    for (int i = 0; i < m_number; ++i) {
        std::cout << "\"" << m_names[i] << "\" ";
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
GModelTemporalRegistry::GModelTemporalRegistry(const GModelTemporalRegistry& registry)
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
GModelTemporalRegistry::~GModelTemporalRegistry(void)
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
GModelTemporalRegistry& GModelTemporalRegistry::operator= (const GModelTemporalRegistry& registry)
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
 * @brief Allocate temporal model of given name
 *
 * @param[in] name Temporal model name.
 * @return Pointer to temporal model (NULL if name is not registered).
 *
 * Returns a pointer to a temporal model instance of the specified name.
 * If the name has not been found in the registry, a NULL pointer is
 * returned.
 ***************************************************************************/
GModelTemporal* GModelTemporalRegistry::alloc(const std::string& name) const
{
    // Initialise temporal model
    GModelTemporal* model = NULL;

    // Search for model in registry
    for (int i = 0; i < m_number; ++i) {
        if (m_names[i] == name) {
            model = m_models[i]->clone();
            break;
        }
    }    

    // Return temporal model
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
std::string GModelTemporalRegistry::name(const int& index) const
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_NAME, index, 0, size()-1);
    }
    #endif

    // Return name
    return (m_names[index]);
}


/***********************************************************************//**
 * @brief Print registry information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return Registry content.
 ***************************************************************************/
std::string GModelTemporalRegistry::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelTemporalRegistry ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of models"));
        result.append(gammalib::str(m_number));

        // Append models
        for (int i = 0; i < m_number; ++i) {
            result.append("\n"+gammalib::parformat(m_names[i]));
            result.append(m_models[i]->type());
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
void GModelTemporalRegistry::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] registry Registry.
 ***************************************************************************/
void GModelTemporalRegistry::copy_members(const GModelTemporalRegistry& registry)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelTemporalRegistry::free_members(void)
{
    // Return
    return;
}
