/***************************************************************************
 *      GModelSpectralRegistry.cpp  -  Spectral model registry class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelSpectralRegistry.cpp
 * @brief GModelSpectralRegistry class interface definition
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GModelSpectralRegistry.hpp"
#include "GTools.hpp"

/* __ Static members _____________________________________________________ */
int                    GModelSpectralRegistry::m_number(0);
std::string*           GModelSpectralRegistry::m_names(0);
const GModelSpectral** GModelSpectralRegistry::m_models(0);

/* __ Method name definitions ____________________________________________ */
#define G_NAME                           "GModelSpectralRegistry::name(int&)"

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
GModelSpectralRegistry::GModelSpectralRegistry(void)
{
    // Initialise private members for clean destruction
    init_members();

    // Debug option: Show actual registry
    #if G_DEBUG_REGISTRY
    std::cout << "GModelSpectralRegistry(void): ";
    for (int i = 0; i < m_number; ++i)
        std::cout << "\"" << m_names[i] << "\" ";
    std::cout << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Model constructor
 ***************************************************************************/
GModelSpectralRegistry::GModelSpectralRegistry(const GModelSpectral* model)
{
    // Initialise private members for clean destruction
    init_members();

    // Debug option: Notify new registry
    #if G_DEBUG_REGISTRY
    std::cout << "GModelSpectralRegistry(const GModelSpatial*): ";
    std::cout << "add \"" << model->type() << "\" to registry." << std::endl;
    #endif
    
    // Allocate new old registry
    std::string*           new_names  = new std::string[m_number+1];
    const GModelSpectral** new_models = new const GModelSpectral*[m_number+1];

    // Save old registry
    for (int i = 0; i < m_number; ++i) {
        new_names[i]  = m_names[i];
        new_models[i] = m_models[i];
    }

    // Add new model to registry
    new_names[m_number]  = model->type();
    new_models[m_number] = model;

    // Delete old registry
    if (m_names  != NULL) delete [] m_names;
    if (m_models != NULL) delete [] m_models;

    // Set pointers on new registry
    m_names  = new_names;
    m_models = new_models;

    // Increment number of models in registry
    m_number++;

    // Debug option: Show actual registry
    #if G_DEBUG_REGISTRY
    std::cout << "GModelSpectralRegistry(const GModelSpatial*): ";
    for (int i = 0; i < m_number; ++i)
        std::cout << "\"" << m_names[i] << "\" ";
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
GModelSpectralRegistry::GModelSpectralRegistry(const GModelSpectralRegistry& registry)
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
GModelSpectralRegistry::~GModelSpectralRegistry(void)
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
 ***************************************************************************/
GModelSpectralRegistry& GModelSpectralRegistry::operator= (const GModelSpectralRegistry& registry)
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
 * @brief Allocate spectral model of given type
 *
 * @param[in] type Spectral model type.
 *
 * Returns a pointer to a void spectral model instance of the specified
 * type. If the type has not been found in the registry, a NULL pointer is
 * returned.
 ***************************************************************************/
GModelSpectral* GModelSpectralRegistry::alloc(const std::string& type) const
{
    // Initialise spectral model
    GModelSpectral* model = NULL;

    // Search for model in registry
    for (int i = 0; i < m_number; ++i) {
        if (m_names[i] == type) {
            model = m_models[i]->clone();
            break;
        }
    }    

    // Return spectral model
    return model;
}


/***********************************************************************//**
 * @brief Returns model name
 *
 * @param[in] index Model index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Model index is out of range.
 ***************************************************************************/
std::string GModelSpectralRegistry::name(const int& index) const
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size())
        throw GException::out_of_range(G_NAME, index, 0, size()-1);
    #endif

    // Return name
    return (m_names[index]);
}


/***********************************************************************//**
 * @brief Print registry information
 ***************************************************************************/
std::string GModelSpectralRegistry::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GModelSpectralRegistry ===");
    result.append("\n"+parformat("Number of models")+str(m_number));

    // Append models
    for (int i = 0; i < m_number; ++i) {
        result.append("\n"+parformat(m_names[i]));
        result.append(m_models[i]->type());
    }

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
void GModelSpectralRegistry::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] registry Registry.
 ***************************************************************************/
void GModelSpectralRegistry::copy_members(const GModelSpectralRegistry& registry)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralRegistry::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                  Friends                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] registry Spectral model registry.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GModelSpectralRegistry& registry)
{
     // Write registry in output stream
    os << registry.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] registry Spectral model registry.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GModelSpectralRegistry& registry)
{
    // Write registry into logger
    log << registry.print();

    // Return logger
    return log;
}
