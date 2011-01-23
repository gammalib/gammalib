/***************************************************************************
 *                  GModels.cpp  -  Model container class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModels.cpp
 * @brief GModels class implementation.
 * @author J. Knodlseder
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
#define G_ACCESS                                  "GModels::operator() (int)"
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
GModels::GModels(void) : GOptimizerPars()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] models Models.
 ***************************************************************************/
GModels::GModels(const GModels& models) : GOptimizerPars(models)
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
 * @param[in] filename XML file.
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
 * @brief Returns pointer to model
 *
 * @param[in] index Model index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Model index is out of range.
 ***************************************************************************/
GModel* GModels::operator() (int index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size())
        throw GException::out_of_range(G_ACCESS, index, 0, size()-1);
    #endif

    // Return pointer
    return m_models[index];
}


/***********************************************************************//**
 * @brief Returns pointer to model
 *
 * @param[in] index Model index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Model index is out of range.
 ***************************************************************************/
const GModel* GModels::operator() (int index) const
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size())
        throw GException::out_of_range(G_ACCESS, index, 0, size()-1);
    #endif

    // Return pointer
    return (m_models[index]);
}


/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] models Models which should be assigned.
 ***************************************************************************/
GModels& GModels::operator= (const GModels& models)
{ 
    // Execute only if object is not identical
    if (this != &models) {

        // Copy base class members
        this->GOptimizerPars::operator=(models);

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


/*==========================================================================
 =                                                                         =
 =                            Public methods                               =
 =                                                                         =
 ==========================================================================*/

 /***********************************************************************//**
 * @brief Clear object.
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GModels::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GOptimizerPars::free_members();

    // Initialise members
    this->GOptimizerPars::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 ***************************************************************************/
GModels* GModels::clone(void) const
{
    return new GModels(*this);
}


/***********************************************************************//**
 * @brief Append model to container
 *
 * @param[in] model Model to be added.
 ***************************************************************************/
void GModels::append(const GModel& model)
{
    // Append model
    m_models.push_back(model.clone());

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load models from XML file.
 *
 * @param[in] file Name of XML file.
 ***************************************************************************/
void GModels::load(const std::string& filename)
{
    // Load XML document
    GXml xml(filename);

    // Read models from XML document
    read(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save models into XML file.
 *
 * @param[in] file Name of XML file.
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
 * Read models from the first source library found in the XML document. It is
 * assumed that each model is composed of a spectral and a spatial model
 * (Fermi-LAT style). The decoding of the spatial and spectral XML elements
 * is done within a GModel constructor.
 *
 * @todo Sources names are not verified so far for uniqueness. This would be
 *       required to achieve an unambiguous update of parameters in an already
 *       existing XML file when using the write method. Also, no control is
 *       performed that only a single spectral and spatial component exists in
 *       each source definition.
 ***************************************************************************/
void GModels::read(const GXml& xml)
{
    // Get pointer on source library
    GXmlElement* lib = xml.element("source_library", 0);

    // Loop over all sources
    int n = lib->elements("source");
    for (int i = 0; i < n; ++i) {

        // Get pointer on source
        GXmlElement* src = (GXmlElement*)lib->element("source", i);

        // Get model type
        std::string type = src->attribute("type");

        // Get model
        GModelRegistry registry;
        GModel*        ptr = registry.alloc(type);

        // If model if valid then read model from XML file
        if (ptr != NULL)
            ptr->read(*src);

        // ... otherwise throw an exception
        else
            throw GException::model_invalid(G_READ, type);

        // Append model
        append(*ptr);

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
        xml.append(new GXmlElement("source_library title=\"source library\""));
    }
    
    // Get pointer on source library
    GXmlElement* lib = xml.element("source_library", 0);

    // Write all sources into library
    for (int i = 0; i < size(); ++i)
        m_models[i]->write(*lib);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns value of source model
 *
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 *
 * This method evaluates the factorized source model at a given set of
 * parameters.
 ***************************************************************************/
/*
double GModels::value(const GSkyDir& srcDir, const GEnergy& srcEng,
                      const GTime& srcTime)
{
    // Initialise value
    double value = 0.0;

    // Evaluate function for all models
    for (int i = 0; i < m_elements; ++i)
        value += m_model[i].value(srcDir, srcEng, srcTime);

    // Return
    return value;
}
*/


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 ***************************************************************************/
double GModels::eval(const GEvent& event, const GObservation& obs)
{
    // Initialise function value
    double value = 0.0;

    // Evaluate function for all models
    for (int i = 0; i < size(); ++i)
        value += m_models[i]->eval(event, obs);

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 ***************************************************************************/
double GModels::eval_gradients(const GEvent& event, const GObservation& obs)
{
    // Initialise function value
    double value = 0.0;

    // Evaluate function for all models
    for (int i = 0; i < size(); ++i)
        value += m_models[i]->eval_gradients(event, obs);

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Print models
 ***************************************************************************/
std::string GModels::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GModels ===");
    result.append("\n"+parformat("Number of models")+str(size()));
    result.append("\n"+parformat("Number of parameters")+str(m_npars));

    // Append models
    for (int i = 0; i < size(); ++i)
        result.append("\n"+m_models[i]->print());

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
 * @param[in] models Models.
 ***************************************************************************/
void GModels::copy_members(const GModels& models)
{
    // Copy members
    m_models = models.m_models;

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModels::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set parameter pointers
 *
 * Gathers all parameter pointers from the models into a linear array of
 * GModelPar pointers. This allows the optimizer to see all parameters in a
 * linear array.
 *
 * @todo Replace GModelPar** by std::vector<GModelPar*> in GOptimizerPars
 ***************************************************************************/
void GModels::set_pointers(void)
{
    // Delete old pointers (if they exist)
    if (m_par != NULL) delete [] m_par;
    m_par = NULL;

    // Determine total number of parameters
    m_npars = 0;
    for (int i = 0; i < size(); ++i)
        m_npars += m_models[i]->size();

    // Continue only if there are parameters
    if (m_npars > 0) {

        // Allocate parameter pointers
        m_par = new GModelPar*[m_npars];

        // Initialise pointer on pointer array
        GModelPar** ptr = m_par;

        // Gather all pointers
        for (int i = 0; i < size(); ++i) {
            for (int k = 0; k < m_models[i]->size(); ++k)
                *ptr++ = &((*m_models[i])(k));
        }

    } // endif: there were model parameters

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
 * @param[in] models Models.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GModels& models)
{
     // Write models in output stream
    os << models.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] models Models.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GModels& models)
{
    // Write models into logger
    log << models.print();

    // Return logger
    return log;
}
