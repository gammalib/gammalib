/***************************************************************************
 *                      GModel.cpp  -  Model class                         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModel.cpp
 * @brief GModel class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GModel.hpp"
#include "GModelTemporalConst.hpp"
#include "GModelSpatialPtsrc.hpp"
#include "GModelSpatialConst.hpp"
#include "GModelSpatialCube.hpp"
#include "GModelSpectralPlaw.hpp"
#include "GModelSpectralPlaw2.hpp"
#include "GModelSpectralExpPlaw.hpp"
#include "GModelSpectralFunc.hpp"
#include "GModelSpectralConst.hpp"
#include "GEvent.hpp"
#include "GPointing.hpp"
#include "GResponse.hpp"
#include "GObservation.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ACCESS                                   "GModel::operator() (int)"
#define G_XML_SPATIAL                     "GModel::xml_spatial(GXmlElement&)"
#define G_XML_SPECTRAL                   "GModel::xml_spectral(GXmlElement&)"
#define G_SPATIAL                "GModel::spatial(GEvent&, GEnergy&, GTime&,"\
                                                       " GObservation, bool)"
#define G_SPECTRAL   "GModel::spectral(GEvent&, GTime&, GObservation&, bool)"
#define G_TEMPORAL           "GModel::temporal(GEvent&, GObservation&, bool)"

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
GModel::GModel(void)
{
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Model from which the instance should be built.
 ***************************************************************************/
GModel::GModel(const GModel& model)
{
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(model);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Construct from spatial and spectral models
 *
 * @param[in] spatial Spatial model.
 * @param[in] spectral Spectral model.
 ***************************************************************************/
GModel::GModel(const GModelSpatial& spatial, const GModelSpectral& spectral)
{
    // Initialise private members for clean destruction
    init_members();

    // Allocate constant
    GModelTemporalConst temporal;

    // Clone spatial and spectral models
    m_spatial  = spatial.clone();
    m_spectral = spectral.clone();
    m_temporal = temporal.clone();

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Construct from spatial and spectral XML elements
 *
 * @param[in] spatial Spatial XML element.
 * @param[in] spectral Spectral XML element.
 ***************************************************************************/
GModel::GModel(const GXmlElement& spatial, const GXmlElement& spectral)
{
    // Initialise private members for clean destruction
    init_members();

    // Allocate constant
    GModelTemporalConst temporal;

    // Clone spatial and spectral models
    m_spatial  = xml_spatial(spatial);
    m_spectral = xml_spectral(spectral);
    m_temporal = temporal.clone();

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GModel::~GModel(void)
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
 * @brief Returns reference to model parameter
 *
 * @param[in] index Parameter index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Parameter index is out of range.
 ***************************************************************************/
GModelPar& GModel::operator() (int index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= m_npars)
        throw GException::out_of_range(G_ACCESS, index, 0, m_npars-1);
    #endif

    // Return pointer
    return *(m_par[index]);
}


/***********************************************************************//**
 * @brief Returns reference to model parameter (const version)
 *
 * @param[in] index Parameter index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Parameter index is out of range.
 ***************************************************************************/
const GModelPar& GModel::operator() (int index) const
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= m_npars)
        throw GException::out_of_range(G_ACCESS, index, 0, m_npars-1);
    #endif

    // Return pointer
    return *(m_par[index]);
}


/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] model Model which should be assigned.
 ***************************************************************************/
GModel& GModel::operator= (const GModel& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(model);

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
 * @brief Clear instance
 *
 * This method properly resets the instance to an initial state.
 ***************************************************************************/
void GModel::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 ***************************************************************************/
GModel* GModel::clone(void) const
{
    return new GModel(*this);
}


/***********************************************************************//**
 * @brief Set instruments to which model applies
 *
 * @param[in] instruments String of instruments.
 *
 * Sets the instruments to which the model applies from a comma separated
 * list of strings. If the instrument string is empty the model is considered
 * to apply to all instruments.
 ***************************************************************************/
void GModel::instruments(const std::string& instruments)
{
    // Clear instruments vector
    m_instruments.clear();

    // Extract instruments
    std::vector<std::string> inst = split(instruments, ",");

    // Attach all instruments
    for (int i = 0; i < inst.size(); ++i)
        m_instruments.push_back(toupper(strip_whitespace(inst[i])));

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
 * This method evaluates the factorized source model for a given true photon
 * arrival direction, true photon energy and true photon arrival time.
 ***************************************************************************/
double GModel::value(const GSkyDir& srcDir, const GEnergy& srcEng,
                     const GTime& srcTime)
{
    // Initialise source model
    double source = 1.0;

    // Evaluate source model
    if (m_spatial  != NULL) source *= m_spatial->eval(srcDir);
    if (m_spectral != NULL) source *= m_spectral->eval(srcEng);
    if (m_temporal != NULL) source *= m_temporal->eval(srcTime);

    // Return
    return source;
}


/***********************************************************************//**
 * @brief Returns parameter gradients of source model
 *
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 *
 * This method returns the parameter gradients of the factorized source model
 * for a given true photon arrival direction, true photon energy and true
 * photon arrival time.
 ***************************************************************************/
GVector GModel::gradients(const GSkyDir& srcDir, const GEnergy& srcEng,
                          const GTime& srcTime)
{
    // Evaluate source model gradients
    if (m_spatial  != NULL) m_spatial->eval_gradients(srcDir);
    if (m_spectral != NULL) m_spectral->eval_gradients(srcEng);
    if (m_temporal != NULL) m_temporal->eval_gradients(srcTime);

    // Set vector of gradients
    GVector gradients;
    if (size() > 0) {
        gradients = GVector(size());
        for (int i = 0; i < size(); ++i)
            gradients(i) = m_par[i]->gradient();
    }

    // Return gradients
    return gradients;
}


/***********************************************************************//**
 * @brief Evaluate model
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 *
 * This method evaluates the source model for a given event within a given
 * observation.
 ***************************************************************************/
double GModel::eval(const GEvent& event, const GObservation& obs)
{
    // Evaluate function
    double value = temporal(event, obs, false);

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate model and parameter gradients
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 *
 * This method evaluates the source model and model gradients for a given
 * event within a given observation.
 ***************************************************************************/
double GModel::eval_gradients(const GEvent& event, const GObservation& obs)
{
    // Evaluate function
    double value = temporal(event, obs, true);

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Verifies if model is valid for a given instrument
 *
 * @param[in] name Instrument name.
 *
 * Checks if specified instrument name is in list of applicable instruments.
 * If the list of applicable instruments is empty the model applies to all
 * possible instruments.
 ***************************************************************************/
bool GModel::isvalid(const std::string& name) const
{
    // Initialise validity
    bool valid = true;

    // Check if model applies to instrument
    if (m_instruments.size() > 0) {

        // Convert instrument name to upper case
        std::string uname = toupper(name);

        // Initialise validity flag
        valid = false;

        // Check if name is in list
        for (int i = 0; i < m_instruments.size(); ++i) {
            if (uname == m_instruments[i]) {
                valid = true;
                break;
            }
        }

    }

    // Return
    return valid;
}


/***********************************************************************//**
 * @brief Print model information
 ***************************************************************************/
std::string GModel::print(void) const
{
    // Initialise result string
    std::string result;

    // Determine number of parameters per type
   int n_spatial  = (m_spatial  != NULL) ? m_spatial->size() : 0;
   int n_spectral = (m_spectral != NULL) ? m_spectral->size() : 0;
   int n_temporal = (m_temporal != NULL) ? m_temporal->size() : 0;

    // Append header
    result.append("=== GModel ===");
    result.append("\n"+parformat("Name")+name());
    result.append("\n"+parformat("Instruments"));
    if (m_instruments.size() > 0) {
        for (int i = 0; i < m_instruments.size(); ++i) {
            if (i > 0)
                result.append(", ");
            result.append(m_instruments[i]);
        }
    }
    else
        result.append("all");
    result.append("\n"+parformat("Number of parameters")+str(size()));
    result.append("\n"+parformat("Number of spatial par's")+str(n_spatial));
    int i;
    for (i = 0; i < n_spatial; ++i) {
        result.append("\n"+(*this)(i).print());
    }
    result.append("\n"+parformat("Number of spectral par's")+str(n_spectral));
    for (; i < n_spatial+n_spectral; ++i) {
        result.append("\n"+(*this)(i).print());
    }
    result.append("\n"+parformat("Number of temporal par's")+str(n_temporal));
    for (; i < n_spatial+n_spectral+n_temporal; ++i) {
        result.append("\n"+(*this)(i).print());
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
void GModel::init_members(void)
{
    // Initialise members
    m_name.clear();
    m_instruments.clear();
    m_npars    = 0;
    m_par      = NULL;
    m_spatial  = NULL;
    m_spectral = NULL;
    m_temporal = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model GModel members which should be copied.
 ***************************************************************************/
void GModel::copy_members(const GModel& model)
{
    // Copy attributes
    m_name        = model.m_name;
    m_instruments = model.m_instruments;
    m_npars       = model.m_npars;

    // Clone spatial and spectral models
    m_spatial  = (model.m_spatial  != NULL) ? model.m_spatial->clone()  : NULL;
    m_spectral = (model.m_spectral != NULL) ? model.m_spectral->clone() : NULL;
    m_temporal = (model.m_temporal != NULL) ? model.m_temporal->clone() : NULL;

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModel::free_members(void)
{
    // Free memory
    if (m_par      != NULL) delete [] m_par;
    if (m_spatial  != NULL) delete m_spatial;
    if (m_spectral != NULL) delete m_spectral;
    if (m_temporal != NULL) delete m_temporal;

    // Signal free pointers
    m_par      = NULL;
    m_spatial  = NULL;
    m_spectral = NULL;
    m_temporal = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set parameter pointers
 *
 * Gathers all parameter pointers from the model.
 ***************************************************************************/
void GModel::set_pointers(void)
{
    // Delete old pointers (if they exist)
    if (m_par != NULL) delete [] m_par;
    m_par = NULL;

    // Determine the number of parameters
    int n_spatial  = (m_spatial  != NULL) ? m_spatial->size() : 0;
    int n_spectral = (m_spectral != NULL) ? m_spectral->size() : 0;
    int n_temporal = (m_temporal != NULL) ? m_temporal->size() : 0;
    m_npars        = n_spatial + n_spectral + n_temporal;

    // Continue only if there are parameters
    if (m_npars > 0) {

        // Allocate parameter pointers
        m_par = new GModelPar*[m_npars];

        // Initialise pointer on pointer array
        GModelPar** ptr = m_par;

        // Gather spatial parameter pointers
        for (int i = 0; i < n_spatial; ++i)
            *ptr++ = &((*m_spatial)(i));

        // Gather spectral parameters
        for (int i = 0; i < n_spectral; ++i)
            *ptr++ = &((*m_spectral)(i));

        // Gather temporal parameters
        for (int i = 0; i < n_temporal; ++i)
            *ptr++ = &((*m_temporal)(i));

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Construct spatial model from XML element
 *
 * @param[in] spatial XML element containing spatial model information.
 *
 * @exception GException::model_invalid_spatial
 *            Invalid spatial model type encountered.
 ***************************************************************************/
GModelSpatial* GModel::xml_spatial(const GXmlElement& spatial) const
{
    // Initialize spatial pointer
    GModelSpatial* ptr = NULL;

    // Get spatial model type
    std::string type = spatial.attribute("type");

    // Type=SkyDirFunction
    if (type == "SkyDirFunction")
        ptr = new GModelSpatialPtsrc(spatial);

    // Type=SpatialMap
    else if (type == "SpatialMap")
        throw GException::model_invalid_spatial(G_XML_SPATIAL, type);

    // Type=MapCubeFunction
    else if (type == "MapCubeFunction")
        ptr = new GModelSpatialCube(spatial);

    // Type=ConstantValue
    else if (type == "ConstantValue")
        ptr = new GModelSpatialConst(spatial);

    // Unknown spatial model type
    else
        throw GException::model_invalid_spatial(G_XML_SPATIAL, type);

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Construct spectral model from XML element
 *
 * @param[in] spectral XML element containing spectral model information.
 *
 * @exception GException::model_invalid_spectral
 *            Invalid spatial model type encountered.
 ***************************************************************************/
GModelSpectral* GModel::xml_spectral(const GXmlElement& spectral) const
{
    // Initialize spectral pointer
    GModelSpectral* ptr = NULL;

    // Get spatial model type
    std::string type = spectral.attribute("type");

    // Type=PowerLaw
    if (type == "PowerLaw")
        ptr = new GModelSpectralPlaw(spectral);

    // Type=PowerLaw2
    else if (type == "PowerLaw2")
        ptr = new GModelSpectralPlaw2(spectral);

    // Type=ExpCutoff
    else if (type == "ExpCutoff")
        ptr = new GModelSpectralExpPlaw(spectral);

    // Type=FileFunction
    else if (type == "FileFunction")
        ptr = new GModelSpectralFunc(spectral);

    // Type=ConstantValue
    else if (type == "ConstantValue")
        ptr = new GModelSpectralConst(spectral);

    // Unknown spectral model type
    else
        throw GException::model_invalid_spectral(G_XML_SPECTRAL, type);

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Returns spatial model component
 *
 * @param[in] event Observed event.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observation.
 * @parem[in] grad Compute also model gradients? (default=false)
 *
 * @exception GException::no_response
 *            Observation has no valid instrument response
 *
 * This method computes the spatial model component for a given true photon
 * energy and arrival time.
 ***************************************************************************/
double GModel::spatial(const GEvent& event,
                       const GEnergy& srcEng, const GTime& srcTime,
                       const GObservation& obs, bool grad)
{
    // Initialise result
    double value = 0.0;

    // Continue only if the model has a spatial component
    if (m_spatial != NULL) {

        // Get response function
        GResponse* rsp = obs.response();
        if (rsp == NULL)
            throw GException::no_response(G_SPATIAL);

        // Get IRF value
        double irf = rsp->irf(event, *this, srcEng, srcTime, obs);

        // Compute source model
        double source = 1.0;

        // Evaluate model in all components
        if (grad) {

            // Evaluate source model
            if (m_spectral != NULL) source *= m_spectral->eval_gradients(srcEng);
            if (m_temporal != NULL) source *= m_temporal->eval_gradients(srcTime);

            // Set value
            value = source * irf;

            // Set gradients
            for (int i = 0; i < m_npars; ++i)
                m_par[i]->gradient(m_par[i]->gradient() * irf);

        }
        else {

            // Evaluate source model
            if (m_spectral != NULL) source *= m_spectral->eval(srcEng);
            if (m_temporal != NULL) source *= m_temporal->eval(srcTime);

            // Set value
            value = source * irf;
        }

    } // endif: Gamma-ray source model had a spatial component

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Perform integration over spectral component
 *
 * @param[in] event Observed event.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observation.
 * @parem[in] grad Compute also model gradients (default=false).
 *
 * @exception GException::no_response
 *            Observation has no valid instrument response
 *
 * This method integrates the source model over the spectral component. If
 * the response function has no energy dispersion then no spectral
 * integration is needed and the observed photon energy is identical to the
 * true photon energy.
 *
 * @todo Needs implementation of spectral integration to handle energy
 *       dispersion.
 ***************************************************************************/
double GModel::spectral(const GEvent& event, const GTime& srcTime,
                        const GObservation& obs, bool grad)
{
    // Initialise result
    double value = 0.0;

    // Get response function
    GResponse* rsp = obs.response();
    if (rsp == NULL)
        throw GException::no_response(G_SPECTRAL);

    // Determine if energy integration is needed
    bool integrate = rsp->hasedisp();

    // Case A: Integraion
    if (integrate) {
        std::cout << "GModel::spectral: Integration not implemented." << std::endl;
    }

    // Case B: No integration (assume no energy dispersion)
    else
        value = spatial(event, event.energy(), srcTime, obs, grad);

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Perform integration over temporal component
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 * @parem[in] grad Compute also model gradients (default=false).
 *
 * @exception GException::no_response
 *            Observation has no valid instrument response
 *
 * This method integrates the source model over the temporal component. If
 * the response function has no time dispersion then no temporal integration
 * is needed and the observed photon arrival time is identical to the true
 * photon arrival time.
 *
 * @todo Needs implementation of temporal integration to handle time
 *       dispersion.
 ***************************************************************************/
double GModel::temporal(const GEvent& event, const GObservation& obs, bool grad)
{
    // Initialise result
    double value = 0.0;

    // Get response function
    GResponse* rsp = obs.response();
    if (rsp == NULL)
        throw GException::no_response(G_TEMPORAL);

    // Determine if time integration is needed
    bool integrate = rsp->hastdisp();

    // Case A: Integraion
    if (integrate) {
        std::cout << "GModel::temporal: Integration not implemented." << std::endl;
    }

    // Case B: No integration (assume no time dispersion)
    else
        value = spectral(event, event.time(), obs, grad);

    // Return value
    return value;
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
 * @param[in] model Model.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GModel& model)
{
     // Write spectrum in output stream
    os << model.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] model Model.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GModel& model)
{
    // Write spectrum into logger
    log << model.print();

    // Return logger
    return log;
}
