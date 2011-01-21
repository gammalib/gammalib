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
#include "GXml.hpp"
#include "GXmlElement.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ACCESS                                  "GModels::operator() (int)"

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
    // Initialise private members for clean destruction
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
    // Initialise private members for clean destruction
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
    // Initialise private members for clean destruction
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
    if (index < 0 || index >= m_elements)
        throw GException::out_of_range(G_ACCESS, index, 0, m_elements-1);
    #endif

    // Return reference
    return &(m_model[index]);
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
    if (index < 0 || index >= m_elements)
        throw GException::out_of_range(G_ACCESS, index, 0, m_elements-1);
    #endif

    // Return pointer
    return &(m_model[index]);
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

        // Initialise private members for clean destruction
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
    // Allocate fresh memory for models
    GModel* new_model = new GModel[m_elements+1];

    // If we have already models then copy them over to the array
    if (m_elements > 0) {
        for (int i = 0; i < m_elements; ++i)
            new_model[i] = m_model[i];
    }

    // Copy the model in the last elements of the array
    new_model[m_elements] = model;

    // Release old model array
    if (m_model != NULL) delete [] m_model;

    // Attach new model array to this object
    m_model = new_model;

    // Increment number of models
    m_elements++;

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
 * Read models from the first source library found in the XML document. It is
 * assumed that each model is composed of a spectral and a spatial model
 * (Fermi-LAT style). The decoding of the spatial and spectral XML elements
 * is done within a GModel constructor.
 *
 * @todo Sources names are not verified so far for uniqueness. This would be
 * required to achieve an unambiguous update of parameters in an already
 * existing XML file when using the write method. Also, no control is
 * performed that only a single spectral and spatial component exists in
 * each source definition.
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

        // Get pointers on spectrum and spatial model
        GXmlElement* spec = (GXmlElement*)src->element("spectrum", 0);
        GXmlElement* spat = (GXmlElement*)src->element("spatialModel", 0);

        // Build model from GXml elements
        GModel model(*spat, *spec);

        // Set model name
        model.name(src->attribute("name"));

        // Set instruments
        model.instruments(src->attribute("instrument"));

        // Append model
        append(model);

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
 *
 * @todo No detailed check is performed to verify the coherence of any
 * existing source elements.
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
    for (int i = 0; i < m_elements; ++i) {

        // Initialise pointer on source
        GXmlElement* src = NULL;

        // Search corresponding source
        int n = lib->elements("source");
        for (int k = 0; k < n; ++k) {
            GXmlElement* element = (GXmlElement*)lib->element("source", k);
            if (element->attribute("name") == m_model[i].name()) {
                src = element;
                break;
            }
        }

        // If no source with corresponding name was found then append one
        if (src == NULL) {
            src = new GXmlElement("source");
            src->attribute("name") = m_model[i].name();
            src->append(new GXmlElement("spectrum"));
            src->append(new GXmlElement("spatialModel"));
            lib->append(src);
        }

        // Set source attributes
        if (m_model[i].spatial()->type() == "PointSource")
            src->attribute("type", "PointSource");
        else
            src->attribute("type", "DiffuseSource");
        src->attribute("name", m_model[i].name());
        std::string instruments = m_model[i].instruments();
        if (instruments.length() > 0)
            src->attribute("instrument", instruments);

        // Get pointers on spectrum and spatial model
        GXmlElement* spec = (GXmlElement*)src->element("spectrum", 0);
        GXmlElement* spat = (GXmlElement*)src->element("spatialModel", 0);

        // Write spectral and spatial information
        m_model[i].spectral()->write(*spec);
        m_model[i].spatial()->write(*spat);

    } // endfor: wrote all sources into library

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
    for (int i = 0; i < m_elements; ++i)
        value += m_model[i].eval(event, obs);

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

    // Evaluate function and gradients for all models
    for (int i = 0; i < m_elements; ++i)
        value += m_model[i].eval_gradients(event, obs);

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Return simulated list of photons
 *
 * @param[in] area Simulation surface area (cm2).
 * @param[in] dir Centre of simulation cone.
 * @param[in] radius Radius of simulation cone (deg).
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 * @param[in] tmin Minimum photon arrival time.
 * @param[in] tmax Maximum photon arrival time.
 * @param[in] ran Random number generator. 
 *
 * This method returns a list of photons that has been derived by Monte Carlo
 * simulation from the models.
 *
 * @todo Implement time ordering.
 ***************************************************************************/
GPhotons GModels::mc(const double& area,
                     const GSkyDir& dir,  const double&  radius,
                     const GEnergy& emin, const GEnergy& emax,
                     const GTime&   tmin, const GTime&   tmax,
                     GRan& ran)
{
    // Allocate photons
    GPhotons photons;

    // Simulate photons for all models
    for (int i = 0; i < m_elements; ++i) {

        // Get photons for actual model
        GPhotons p = m_model[i].mc(area, dir, radius, emin, emax, tmin, tmax, ran);

        // Reserve new space for photons
        photons.reserve(photons.size() + p.size());
        
        // Add photons to list
        for (int k = 0; k < p.size(); ++k) {
            p[k].mcid(i);
            photons.push_back(p[k]);
        }

    } // endfor: looped over models

    // Return photon list
    return photons;
}


/***********************************************************************//**
 * @brief Print models
 ***************************************************************************/
std::string GModels::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GModels ===\n");
    result.append(parformat("Number of models")+str(size())+"\n");
    result.append(parformat("Number of parameters")+str(m_npars));

    // Append models
    for (int k = 0; k < size(); ++k) {
        const GModel* ptr = (*this)(k);
        result.append("\n"+parformat("Model name")+(*this)(k)->name());
        result.append("\n"+parformat("Model type"));
        result.append(ptr->temporal()->type()+" ");
        result.append(ptr->spectral()->type()+" ");
        result.append(ptr->spatial()->type());
        result.append("\n"+parformat("Instruments"));
        if (ptr->m_instruments.size() > 0) {
            for (int i = 0; i < ptr->m_instruments.size(); ++i) {
                if (i > 0)
                    result.append(", ");
                result.append(ptr->m_instruments[i]);
            }
        }
        else
            result.append("all");
        for (int j = 0; j < ptr->size(); ++j)
            result.append("\n"+(*ptr)(j).print());
    }
    
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
    m_elements = 0;
    m_model    = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] models GModels members which should be copied.
 ***************************************************************************/
void GModels::copy_members(const GModels& models)
{
    // Copy attributes
    m_elements = models.m_elements;

    // If there are models then copy them
    if (m_elements > 0 && models.m_model != NULL) {

        // Allocate models
        m_model = new GModel[m_elements];

        // Copy models
        for (int i = 0; i < m_elements; ++i)
            m_model[i] = models.m_model[i];

    }

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
    // Free memory
    if (m_model != NULL) delete [] m_model;

    // Signal free pointers
    m_model = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set parameter pointers
 *
 * Gathers all parameter pointers from the models into a linear array of
 * GModelPar pointers. This allows the optimizer to see all parameters in a
 * linear array.
 ***************************************************************************/
void GModels::set_pointers(void)
{
    // Delete old pointers (if they exist)
    if (m_par != NULL) delete [] m_par;
    m_par = NULL;

    // Determine total number of parameters
    m_npars = 0;
    for (int k = 0; k < m_elements; ++k)
        m_npars += m_model[k].size();

    // Continue only if there are parameters
    if (m_npars > 0) {

        // Allocate parameter pointers
        m_par = new GModelPar*[m_npars];

        // Initialise pointer on pointer array
        GModelPar** ptr = m_par;

        // Loop over all models
        for (int k = 0; k < m_elements; ++k) {

            // Determine the number of parameters of each type
            int n_spatial  = (m_model[k].spatial()  != NULL) 
                             ? m_model[k].spatial()->size() : 0;
            int n_spectral = (m_model[k].spectral() != NULL) 
                             ? m_model[k].spectral()->size() : 0;
            int n_temporal = (m_model[k].temporal() != NULL) 
                             ? m_model[k].temporal()->size() : 0;

            // Gather spatial parameter pointers
            for (int i = 0; i < n_spatial; ++i)
                *ptr++ = &((*(m_model[k].spatial()))(i));

            // Gather spectral parameters
            for (int i = 0; i < n_spectral; ++i)
                *ptr++ = &((*(m_model[k].spectral()))(i));

            // Gather temporal parameters
            for (int i = 0; i < n_temporal; ++i)
                *ptr++ = &((*(m_model[k].temporal()))(i));

        } // endfor: looped over all models

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
