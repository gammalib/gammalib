/***************************************************************************
 *       GCTAModelSpatialRegistry.cpp - Spatial model registry class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2018 by Juergen Knoedlseder                              *
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
 * @file GCTAModelSpatialRegistry.cpp
 * @brief Spatial model registry class implemenation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GModelPar.hpp"
#include "GCTAModelSpatial.hpp"
#include "GCTAModelSpatialRegistry.hpp"
#include "GException.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ALLOC               "GCTAModelSpatialRegistry::alloc(GXmlElement&)"
#define G_NAME                         "GCTAModelSpatialRegistry::name(int&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_REGISTRY        //!< Dump registry


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTAModelSpatialRegistry::GCTAModelSpatialRegistry(void)
{
    // Initialise private members for clean destruction
    init_members();

    // Debug option: Show actual registry
    #if defined(G_DEBUG_REGISTRY)
    std::cout << "GCTAModelSpatialRegistry(void): ";
    for (int i = 0; i < size(); ++i) {
        std::cout << "\"" << models()[i]->type() << "\" ";
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
GCTAModelSpatialRegistry::GCTAModelSpatialRegistry(const GCTAModelSpatial* model)
{
    // Initialise private members for clean destruction
    init_members();

    // Debug option: Notify new registry
    #if defined(G_DEBUG_REGISTRY)
    std::cout << "GCTAModelSpatialRegistry(const GCTAModelSpatial*): ";
    std::cout << "add \"" << model->type() << "\" to registry." << std::endl;
    #endif

    // Allocate new registry
    const GCTAModelSpatial** new_models = new const GCTAModelSpatial*[size()+1];

    // Save old registry
    for (int i = 0; i < size(); ++i) {
        new_models[i] = models()[i];
    }

    // Add new model to registry
    new_models[size()] = model;

    // Set pointers on new registry
    models().assign(new_models);

    // Increment number of models in registry
    number()++;

    // Debug option: Show actual registry
    #if defined(G_DEBUG_REGISTRY)
    std::cout << "GCTAModelSpatialRegistry(const GCTAModelSpatial*): ";
    for (int i = 0; i < size(); ++i) {
        std::cout << "\"" << models()[i]->type() << "\" ";
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
GCTAModelSpatialRegistry::GCTAModelSpatialRegistry(const GCTAModelSpatialRegistry& registry)
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
GCTAModelSpatialRegistry::~GCTAModelSpatialRegistry(void)
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
GCTAModelSpatialRegistry& GCTAModelSpatialRegistry::operator=(const GCTAModelSpatialRegistry& registry)
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
 * @brief Allocate spatial model that is found in XML element
 *
 * @param[in] xml XML element.
 * @return Pointer to spatial model.
 *
 * @exception GException::invalid_value
 *            No appropriate spatial model found in XML element.
 *
 * Returns a pointer to a spatial model instance that corresponds to the
 * type found in an XML element. If no appropriate model is found the method
 * will throw an exception.
 ***************************************************************************/
GCTAModelSpatial* GCTAModelSpatialRegistry::alloc(const GXmlElement& xml) const
{
    // Initialise spatial model
    GCTAModelSpatial* model = NULL;

    // Search for model type in registry
    for (int i = 0; i < size(); ++i) {
        if (models()[i]->type() == xml.attribute("type")) {
            model = models()[i]->clone();
            break;
        }
    }

    // If no model has been found then throw an exception
    if (model == NULL) {
        std::string msg = "Spatial model of type \""+xml.attribute("type")+
                          "\" not found in registry. Possible spatial model "
                          "types are:";
        for (int i = 0; i < size(); ++i) {
            msg += " \""+models()[i]->type()+"\"";
        }
        msg += ".";
        throw GException::invalid_value(G_ALLOC, msg);
    }

    // Read model from XML element
    model->read(xml);

    // Return spectral model
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
std::string GCTAModelSpatialRegistry::name(const int& index) const
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_NAME, "Spatial model index", index,
                                       size());
    }
    #endif

    // Return name
    return (models()[index]->type());
}


/***********************************************************************//**
 * @brief Print registry information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return Registry content.
 ***************************************************************************/
std::string GCTAModelSpatialRegistry::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAModelSpatialRegistry ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of models"));
        result.append(gammalib::str(size()));

        // NORMAL: Append models
        if (chatter >= NORMAL) {
            for (int i = 0; i < size(); ++i) {
                result.append("\n"+gammalib::parformat(models()[i]->type()));
                for (int k = 0; k < models()[i]->size(); ++k) {
                    if (k > 0) {
                        result.append(", ");
                    }
                    result.append("\"");
                    result.append((*(models()[i]))[k].name());
                    result.append("\"");
                }
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
void GCTAModelSpatialRegistry::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] registry Registry.
 ***************************************************************************/
void GCTAModelSpatialRegistry::copy_members(const GCTAModelSpatialRegistry& registry)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAModelSpatialRegistry::free_members(void)
{
    // Return
    return;
}
