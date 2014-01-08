/***************************************************************************
 *                   GModels.cpp - Model container class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2013 by Juergen Knoedlseder                         *
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
 * @file GModels.cpp
 * @brief Model container class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GModels.hpp"
#include "GModel.hpp"
#include "GModelRegistry.hpp"
#include "GXml.hpp"
#include "GXmlElement.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ACCESS                          "GModels::operator[](std::string&)"
#define G_AT                                              "GModels::at(int&)"
#define G_SET1                                  "GModels::set(int&, GModel&)"
#define G_SET2                          "GModels::set(std::string&, GModel&)"
#define G_APPEND                                   "GModels::append(GModel&)"
#define G_INSERT1                            "GModels::insert(int&, GModel&)"
#define G_INSERT2                    "GModels::insert(std::string&, GModel&)"
#define G_REMOVE1                                     "GModels::remove(int&)"
#define G_REMOVE2                             "GModels::remove(std::string&)"
#define G_EXTEND                                  "GModels::extend(GModels&)"
#define G_READ                                         "GModels::read(GXml&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GModels::GModels(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] models Model container.
 ***************************************************************************/
GModels::GModels(const GModels& models)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(models);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load constructor
 *
 * @param[in] filename XML filename.
 *
 * Constructs model container from an XML file. See the read() method for
 * more information about the expected structure of the XML file.
 ***************************************************************************/
GModels::GModels(const std::string& filename)
{
    // Initialise members
    init_members();

    // Load XML file
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GModels::~GModels(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] models Model container.
 * @return Model container.
 ***************************************************************************/
GModels& GModels::operator=(const GModels& models)
{
    // Execute only if object is not identical
    if (this != &models) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(models);

    } // endif: object was not identical

    // Return
    return *this;
}


/***********************************************************************//**
 * @brief Return pointer to model
 *
 * @param[in] name Model name.
 *
 * @exception GException::model_not_found
 *            Model with specified name not found in container.
 *
 * Returns a pointer to the model with the specified @p name.
 ***************************************************************************/
GModel* GModels::operator[](const std::string& name)
{
    // Get model index
    int index = get_index(name);

    // Throw exception if model name was not found
    if (index == -1) {
        throw GException::model_not_found(G_ACCESS, name);
    }

    // Return pointer
    return m_models[index];
}


/***********************************************************************//**
 * @brief Return pointer to model (const version)
 *
 * @param[in] name Model name.
 *
 * @exception GException::model_not_found
 *            Model with specified name not found in container.
 *
 * Returns a const pointer to the model with the specified @p name.
 ***************************************************************************/
const GModel* GModels::operator[](const std::string& name) const
{
    // Get model index
    int index = get_index(name);

    // Throw exception if model name was not found
    if (index == -1) {
        throw GException::model_not_found(G_ACCESS, name);
    }

    // Return pointer
    return m_models[index];
}


/*==========================================================================
 =                                                                         =
 =                            Public methods                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear object
 *
 * Removes all models from the container.
 ***************************************************************************/
void GModels::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 *
 * @return Pointer to deep copy of model container
 *
 * Makes a deep copy of the model container instance.
 ***************************************************************************/
GModels* GModels::clone(void) const
{
    return new GModels(*this);
}


/***********************************************************************//**
 * @brief Return pointer to model
 *
 * @param[in] index Model index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Model index is out of range.
 *
 * Returns a pointer to the model with the specified @p index.
 ***************************************************************************/
GModel* GModels::at(const int& index)
{
    // Compile option: raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Model index", index, size());
    }
    #endif

    // Return pointer
    return m_models[index];
}


/***********************************************************************//**
 * @brief Return pointer to model (const version)
 *
 * @param[in] index Model index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Model index is out of range.
 *
 * Returns a const pointer to the model with the specified @p index.
 ***************************************************************************/
const GModel* GModels::at(const int& index) const
{
    // Compile option: raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Model index", index, size());
    }
    #endif

    // Return pointer
    return m_models[index];
}


/***********************************************************************//**
 * @brief Set model in container
 *
 * @param[in] index Model index [0,...,size()-1].
 * @param[in] model Model.
 * @return Pointer to deep copy of model.
 *
 * @exception GException::out_of_range
 *            Model index is out of range.
 * @exception GException::invalid_value
 *            Name of model exists already in container.
 *
 * Set model in the container. A deep copy of the model will be made.
 ***************************************************************************/
GModel* GModels::set(const int& index, const GModel& model)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_SET1, index, 0, size()-1);
    }
    #endif

    // Check if a model with specified name does not yet exist
    int inx = get_index(model.name());
    if (inx != -1 && inx != index) {
        std::string msg =
            "Attempt to set model with name \""+model.name()+"\" in model"
            " container at index "+gammalib::str(index)+", but a model with"
            " the same name exists already at index "+gammalib::str(inx)+
            " in the container.\n"
            "Every model in the model container needs a unique name.";
        throw GException::invalid_value(G_SET1, msg);
    }

    // Delete any existing model
    if (m_models[index] != NULL) delete m_models[index];

    // Assign new model by cloning
    m_models[index] = model.clone();

    // Return pointer to model
    return m_models[index];
}


/***********************************************************************//**
 * @brief Set model in container
 *
 * @param[in] name Model name.
 * @param[in] model Model pointer.
 * @return Pointer to deep copy of model.
 *
 * @exception GException::model_not_found
 *            Model with specified name not found in container.
 * @exception GException::invalid_value
 *            Name of model exists already in container.
 *
 * Set model in the container. A deep copy of the model will be made.
 ***************************************************************************/
GModel* GModels::set(const std::string& name, const GModel& model)
{
    // Get parameter index
    int index = get_index(name);

    // Throw exception if parameter name was not found
    if (index == -1) {
        throw GException::model_not_found(G_SET2, name);
    }

    // Check if a model with specified name does not yet exist
    int inx = get_index(model.name());
    if (inx != -1 && inx != index) {
        std::string msg =
            "Attempt to set model with name \""+model.name()+"\" in model"
            " container at index "+gammalib::str(index)+", but a model with"
            " the same name exists already at index "+gammalib::str(inx)+
            " in the container.\n"
            "Every model in the model container needs a unique name.";
        throw GException::invalid_value(G_SET2, msg);
    }

    // Delete any existing model
    if (m_models[index] != NULL) delete m_models[index];

    // Assign new model by cloning
    m_models[index] = model.clone();

    // Return pointer to model
    return m_models[index];
}


/***********************************************************************//**
 * @brief Append model to container
 *
 * @param[in] model Model.
 * @return Pointer to deep copy of model.
 *
 * @exception GException::invalid_value
 *            Name of model exists already in container.
 *
 * Appends model to the container by making a deep copy of the model and
 * storing its pointer.
 ***************************************************************************/
GModel* GModels::append(const GModel& model)
{
    // Check if a model with specified name does not yet exist
    int inx = get_index(model.name());
    if (inx != -1) {
        std::string msg = 
            "Attempt to append model with name \""+model.name()+"\" to model"
            " container, but a model with the same name exists already at"
            " index "+gammalib::str(inx)+" in the container.\n"
            "Every model in the model container needs a unique name.";
        throw GException::invalid_value(G_APPEND, msg);
    }

    // Create deep copy of model
    GModel* ptr = model.clone();

    // Append deep copy of model
    m_models.push_back(ptr);

    // Return pointer to model
    return ptr;
}


/***********************************************************************//**
 * @brief Insert model into container
 *
 * @param[in] index Model index [0,...,size()-1].
 * @param[in] model Model.
 * @return Pointer to deep copy of model.
 *
 * @exception GException::out_of_range
 *            Model index is out of range.
 * @exception GException::invalid_value
 *            Name of model exists already in container.
 *
 * Inserts a @p model into the container before the model with the specified
 * @p index.
 ***************************************************************************/
GModel* GModels::insert(const int& index, const GModel& model)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (is_empty()) {
        if (index > 0) {
            throw GException::out_of_range(G_INSERT1, index, 0, size()-1);
        }
    }
    else {
        if (index < 0 || index >= size()) {
            throw GException::out_of_range(G_INSERT1, index, 0, size()-1);
        }
    }
    #endif

    // Check if a model with specified name does not yet exist
    int inx = get_index(model.name());
    if (inx != -1) {
        std::string msg =
            "Attempt to insert model with name \""+model.name()+"\" in model"
            " container before index "+gammalib::str(index)+", but a model"
            " with the same name exists already at index "+gammalib::str(inx)+
            " in the container.\n"
            "Every model in the model container needs a unique name.";
        throw GException::invalid_value(G_INSERT1, msg);
    }

    // Create deep copy of model
    GModel* ptr = model.clone();

    // Inserts deep copy of model
    m_models.insert(m_models.begin()+index, ptr);

    // Return pointer to model
    return ptr;
}


/***********************************************************************//**
 * @brief Insert model into container
 *
 * @param[in] name Model name.
 * @param[in] model Model.
 * @return Pointer to deep copy of model.
 *
 * @exception GException::model_not_found
 *            Model with specified name not found in container.
 * @exception GException::invalid_value
 *            Name of model exists already in container.
 *
 * Inserts a @p model into the container before the model with the specified
 * @p name.
 ***************************************************************************/
GModel* GModels::insert(const std::string& name, const GModel& model)
{
    // Get parameter index
    int index = get_index(name);

    // Throw exception if parameter name was not found
    if (index == -1) {
        throw GException::model_not_found(G_INSERT2, name);
    }

    // Check if a model with specified name does not yet exist
    int inx = get_index(model.name());
    if (inx != -1) {
        std::string msg =
            "Attempt to insert model with name \""+model.name()+"\" in model"
            " container before index "+gammalib::str(index)+", but a model"
            " with the same name exists already at index "+gammalib::str(inx)+
            " in the container.\n"
            "Every model in the model container needs a unique name.";
        throw GException::invalid_value(G_INSERT2, msg);
    }

    // Create deep copy of model
    GModel* ptr = model.clone();

    // Inserts deep copy of model
    m_models.insert(m_models.begin()+index, ptr);

    // Return pointer to model
    return ptr;
}


/***********************************************************************//**
 * @brief Remove model from container
 *
 * @param[in] index Model index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Model index is out of range.
 *
 * Remove model of specified @p index from container.
 ***************************************************************************/
void GModels::remove(const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_REMOVE1, index, 0, size()-1);
    }
    #endif

    // Delete model
    delete m_models[index];

    // Erase model component from container
    m_models.erase(m_models.begin() + index);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Remove model from container
 *
 * @param[in] name Model name.
 *
 * @exception GException::model_not_found
 *            Model with specified name not found in container.
 *
 * Remove model of specified @p name from container.
 ***************************************************************************/
void GModels::remove(const std::string& name)
{
    // Get parameter index
    int index = get_index(name);

    // Throw exception if parameter name was not found
    if (index == -1) {
        throw GException::model_not_found(G_REMOVE2, name);
    }

    // Delete model
    delete m_models[index];

    // Erase model component from container
    m_models.erase(m_models.begin() + index);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append model container
 *
 * @param[in] models Model container.
 *
 * Append model container to the container.
 ***************************************************************************/
void GModels::extend(const GModels& models)
{
    // Do nothing if model container is empty
    if (!models.is_empty()) {

        // Get size. Note that we extract the size first to avoid an
        // endless loop that arises when a container is appended to
        // itself.
        int num = models.size();

        // Reserve enough space
        reserve(size() + num);

        // Loop over all model components and append pointers to deep copies 
        for (int i = 0; i < num; ++i) {

            // Check if model name does not yet exist
            int inx = get_index(models[i]->name());
            if (inx != -1) {
                std::string msg =
                    "Attempt to append model with name \""+models[i]->name()+
                    "\" to model container, but a model with the same name"
                    " exists already at index "+gammalib::str(inx)+" in the"
                    " container.\n"
                    "Every model in the model container needs a unique name.";
                throw GException::invalid_value(G_EXTEND, msg);
            }

            // Append model to container
            m_models.push_back(models[i]->clone());

        } // endfor: looped over all models

    } // endif: model container was not empty
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Signals if model name exists
 *
 * @param[in] name Model name.
 * @return True if model with specified @p name exists.
 *
 * Searches all model names for a match with the specified @p name. If the
 * specified name has been found, true is returned.
 ***************************************************************************/
bool GModels::contains(const std::string& name) const
{
    // Get model index
    int index = get_index(name);

    // Return
    return (index != -1);
}


/***********************************************************************//**
 * @brief Load models from XML file
 *
 * @param[in] filename XML filename.
 *
 * Loads all models from an XML file. See the read() method for more
 * information about the expected structure of the XML file.
 ***************************************************************************/
void GModels::load(const std::string& filename)
{
    // Clear any existing models
    clear();

    // Load XML document
    GXml xml(filename);

    // Read models from XML document
    read(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save models into XML file
 *
 * @param[in] filename XML filename.
 *
 * Saves all models in the container into an XML file.
 ***************************************************************************/
void GModels::save(const std::string& filename) const
{
    // Declare empty XML document
    GXml xml;

    // Write models into XML file
    write(xml);

    // Save XML document
    xml.save(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read models from XML document
 *
 * @param[in] xml XML document.
 *
 * @exception GException::model_invalid
 *            Invalid model type encountered.
 *
 * Read models from the first source library found in the XML document. The
 * XML document is expected to have the following structure
 *
 *     <source_library title="source library">
 *       <source name="Source1" type="PointSource">
 *         ...
 *       </source>
 *       <source name="Source2" type="DiffuseSource">
 *         ...
 *       </source>
 *     </source_library>
 *
 * Each @p source tag will be interpreted as a model component.
 *
 * @todo Sources names are not verified so far for uniqueness. This would be
 *       required to achieve an unambiguous update of parameters in an already
 *       existing XML file when using the write method.
 ***************************************************************************/
void GModels::read(const GXml& xml)
{
    // Get pointer on source library
    const GXmlElement* lib = xml.element("source_library", 0);

    // Loop over all sources
    int n = lib->elements("source");
    for (int i = 0; i < n; ++i) {

        // Get pointer on source
        const GXmlElement* src = lib->element("source", i);

        // Get model type
        std::string type = src->attribute("type");

        // Get model
        GModelRegistry registry;
        GModel*        ptr = registry.alloc(type);

        // If model if valid then read model from XML file
        if (ptr != NULL) {
            ptr->read(*src);
        }

        // ... otherwise throw an exception
        else {
            throw GException::model_invalid(G_READ, type);
        }

        // Append model
        append(*ptr);

        // Free model (appending clones the model)
        delete ptr;

    } // endfor: looped over all sources

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write models into XML document
 *
 * @param[in] xml XML document.
 *
 * Write models into the first source library that is found in the XML
 * document. In case that no source library exists, one is added to the
 * document.
 ***************************************************************************/
void GModels::write(GXml& xml) const
{
    // If there is no source library then append one
    if (xml.elements("source_library") == 0) {
        xml.append(GXmlElement("source_library title=\"source library\""));
    }

    // Get pointer on source library
    GXmlElement* lib = xml.element("source_library", 0);

    // Write all sources into library
    for (int i = 0; i < size(); ++i) {
        m_models[i]->write(*lib);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return number of model parameters in container
 *
 * @return Number of model parameters in container.
 ***************************************************************************/
int GModels::npars(void) const
{
    // Initialise number of parameters
    int npars = 0;

    // Sum of number of parameters in model
    for (int i = 0; i < size(); ++i) {
        npars += m_models[i]->size();
    }

    // Return
    return npars;
}


/***********************************************************************//**
 * @brief Return optimizer parameter container
 *
 * @return Optimizer parameter container.
 *
 * Returns an optimizer parameter container that is built by extracting
 * all model pointers from the model and attaching them to the parameter
 * container. The optimizer parameter container will thus contains a flat
 * array of a model parameters.
 ***************************************************************************/
GOptimizerPars GModels::pars(void)
{
    // Initialise parameter container
    GOptimizerPars pars;

    // Attach all parameters
    for (int i = 0; i < size(); ++i) {
        GModel* model = m_models[i];
        int     npars = model->size();
        for (int k = 0; k < npars; ++k) {
            pars.attach(&((*model)[k]));
        }
    }

    // Return
    return pars;
}


/***********************************************************************//**
 * @brief Evaluate sum of all models
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 *
 * Evaluates the sum of all models for the specified event and observation.
 ***************************************************************************/
double GModels::eval(const GEvent& event, const GObservation& obs) const
{
    // Initialise function value
    double value = 0.0;

    // Evaluate function for all models
    for (int i = 0; i < size(); ++i) {
        value += m_models[i]->eval(event, obs);
    }

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate sum and gradients of all models
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 *
 * Evaluates the sum and the parameter gradients of all models for the
 * specified event and observation.
 ***************************************************************************/
double GModels::eval_gradients(const GEvent& event,
                               const GObservation& obs) const
{
    // Initialise function value
    double value = 0.0;

    // Evaluate function for all models
    for (int i = 0; i < size(); ++i) {
        value += m_models[i]->eval_gradients(event, obs);
    }

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Print models
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing model container information.
 *
 * Prints all models into a string.
 ***************************************************************************/
std::string GModels::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModels ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of models"));
        result.append(gammalib::str(size()));
        result.append("\n"+gammalib::parformat("Number of parameters"));
        result.append(gammalib::str(npars()));

        // Append models
        for (int i = 0; i < size(); ++i) {
            result.append("\n"+m_models[i]->print(chatter));
        }

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                              Private methods                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModels::init_members(void)
{
    // Initialise members
    m_models.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] models Model container.
 *
 * Makes a copy of all class members. All models are deep copied, and the
 * linear pointer array for parameter access through the GOptimizerPars
 * base class is set.
 ***************************************************************************/
void GModels::copy_members(const GModels& models)
{
    // Copy models
    m_models.clear();
    for (int i = 0; i < models.m_models.size(); ++i) {
        m_models.push_back((models.m_models[i]->clone()));
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 *
 * Deallocates all models. The method loops over the model container and
 * deallocates the memory that has been allocated before.
 ***************************************************************************/
void GModels::free_members(void)
{
    // Free models
    for (int i = 0; i < m_models.size(); ++i) {
        delete m_models[i];
        m_models[i] = NULL;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return model index by name
 *
 * @param[in] name Model name.
 * @return Model index (-1 if model name has not been found)
 *
 * Returns model index based on the specified @p name. If no model with the
 * specified @p name is found the method returns -1.
 ***************************************************************************/
int GModels::get_index(const std::string& name) const
{
    // Initialise index
    int index = -1;

    // Search model with specified name
    for (int i = 0; i < size(); ++i) {
        if (m_models[i]->name() == name) {
            index = i;
            break;
        }
    }

    // Return index
    return index;
}
