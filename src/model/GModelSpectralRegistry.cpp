/***************************************************************************
 *       GModelSpectralRegistry.cpp - Spectral model registry class        *
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
 * @file GModelSpectralRegistry.cpp
 * @brief Spectral model registry class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GModelPar.hpp"
#include "GModelSpectral.hpp"
#include "GModelSpectralRegistry.hpp"
#include "GException.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ALLOC                 "GModelSpectralRegistry::alloc(GXmlElement&)"
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
GModelSpectralRegistry::GModelSpectralRegistry(const GModelSpectral* model)
{
    // Initialise private members for clean destruction
    init_members();

    // Debug option: Notify new registry
    #if G_DEBUG_REGISTRY
    std::cout << "GModelSpectralRegistry(const GModelSpatial*): ";
    std::cout << "add \"" << model->type() << "\" to registry." << std::endl;
    #endif

    // Allocate new registry
    const GModelSpectral** new_models = new const GModelSpectral*[size()+1];

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
    #if G_DEBUG_REGISTRY
    std::cout << "GModelSpectralRegistry(const GModelSpectral*): ";
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
 * @return Reference to registry.
 ***************************************************************************/
GModelSpectralRegistry& GModelSpectralRegistry::operator=(const GModelSpectralRegistry& registry)
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
 * @brief Allocate spectral model that is found in XML element
 *
 * @param[in] xml XML element.
 * @return Pointer to spectral model.
 *
 * @exception GException::invalid_value
 *            No appropriate spectral model found in XML element.
 *
 * Returns a pointer to a spectral model instance that corresponds to the
 * type and parameters found in an XML element. If no appropriate model is
 * found the method will throw an exception.
 ***************************************************************************/
GModelSpectral* GModelSpectralRegistry::alloc(const GXmlElement& xml) const
{
    // Initialise spectral model
    GModelSpectral* model = NULL;

    // Initialise missing parameters string
    std::string missing = "";

    // Search for appropriate model in registry
    for (int i = 0; i < size(); ++i) {

        // If an appropriate model type was found then check the model
        // parameters
        if (models()[i]->type() == xml.attribute("type")) {

            // Check if all required parameters are present. Put all missing
            // parameters as a comma separated list in the "missing" string.
            // If at least one parameter is missing the "found" flag will be
            // set to "false".
            bool found = true;
            for (int k = 0; k < models()[i]->size(); ++k) {
                if (!gammalib::xml_has_par(xml, (*(models()[i]))[k].name())) {
                    if (missing.length() > 0) {
                        missing += " ,";
                    }
                    missing += "\""+(*(models()[i]))[k].name()+"\"";
                    found = false;
                }
            }

            // If parameters are missing then check the next model
            if (!found) {
                continue;
            }

            // No parameters are missing. We thus have the appropriate
            // model and can break now.
            model = models()[i]->clone();
            break;

        } // endif: appropriate type found

    } // endfor: looped over models

    // If no model has been found then throw an exception
    if (model == NULL) {

        // Initialise exception message
        std::string msg = "";

        // If the type has been found then we had missing parameters
        if (missing.length() > 0) {
            msg = "Spectral model of type \""+xml.attribute("type")+ "\" found "
                  "but the following parameters are missing: "+missing;
        }

        // ... otherwise the appropriate type was not found
        else {
            msg = "Spectral model of type \""+xml.attribute("type")+ "\" not "
                  "found in registry. Possible spectral model types are:";
            for (int i = 0; i < size(); ++i) {
                msg += " \""+models()[i]->type()+"\"";
            }
            msg += ".";
        }

        // Throw exception
        throw GException::invalid_value(G_ALLOC, msg);

    } // endif: model was not found

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
std::string GModelSpectralRegistry::name(const int& index) const
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_NAME, "Spectral model index", index,
                                       size());
    }
    #endif

    // Return name
    return (models()[index]->type());
}


/***********************************************************************//**
 * @brief Print registry information
 *
 * @param[in] chatter Chattiness.
 * @return Registry content.
 ***************************************************************************/
std::string GModelSpectralRegistry::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpectralRegistry ===");
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
