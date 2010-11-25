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
#include "GException.hpp"
#include "GModel.hpp"
#include "GModelTemporalConst.hpp"
#include "GModelSpatialPtsrc.hpp"
#include "GModelSpectralPlaw.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_PAR                                              "GModel::par(int)"
#define G_XML_SPATIAL                     "GModel::xml_spatial(GXmlElement&)"
#define G_XML_SPECTRAL                   "GModel::xml_spectral(GXmlElement&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
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
GModel::~GModel()
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
 * @brief Returns pointer to a model parameter
 *
 * @param[in] index Parameter index.
 *
 * @exception GException::out_of_range
 *            Parameter index is out of range.
 ***************************************************************************/
GModelPar* GModel::par(int index) const
{
    // If index is outside boundary then throw an error
    if (index < 0 || index >= m_npars)
        throw GException::out_of_range(G_PAR, index, 0, m_npars-1);

    // Return parameter pointer
    return m_par[index];
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
 * at a given set of parameters.
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
    if (npars() > 0) {
        gradients = GVector(npars());
        for (int i = 0; i < npars(); ++i)
            gradients(i) = m_par[i]->gradient();
    }

    // Return gradients
    return gradients;
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] obsDir Observed photon direction.
 * @param[in] obsEng Observed photon energy.
 * @param[in] obsTime Observed photon arrival time.
 * @param[in] rsp Instrument response function.
 * @param[in] pnt Instrument pointing direction.
 *
 * This method evaluates the source model.
 ***************************************************************************/
double GModel::eval(const GInstDir& obsDir, const GEnergy& obsEng,
                    const GTime& obsTime, const GResponse& rsp,
                    const GPointing& pnt)
{
    // Evaluate function
    double value = temporal(obsDir, obsEng, obsTime, rsp, pnt, false);

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] obsDir Observed photon direction.
 * @param[in] obsEng Observed photon energy.
 * @param[in] obsTime Observed photon arrival time.
 * @param[in] rsp Instrument response function.
 * @param[in] pnt Instrument pointing direction.
 *
 * This method evaluates the source model and sets the model gradients.
 ***************************************************************************/
double GModel::eval_gradients(const GInstDir& obsDir, const GEnergy& obsEng,
                              const GTime& obsTime, const GResponse& rsp,
                              const GPointing& pnt)
{
    // Evaluate function
    double value = temporal(obsDir, obsEng, obsTime, rsp, pnt, true);

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Verifies if model is valid for a given instrument
 *
 * @param[in] name Instrument name for which validity is to be checked.
 *
 * @todo Implement method (dummy method always returns true)
 ***************************************************************************/
bool GModel::isvalid(const std::string& name) const
{
    // Return
    return true;
}


/*==========================================================================
 =                                                                         =
 =                          GModel private methods                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModel::init_members(void)
{
    // Initialise members
    m_name.clear();
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
    m_name  = model.m_name;
    m_npars = model.m_npars;

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
 ***************************************************************************/
void GModel::set_pointers(void)
{
    // Delete old pointers (if they exist)
    if (m_par != NULL) delete [] m_par;
    m_par = NULL;

    // Determine the number of parameters
    int n_spatial  = (m_spatial  != NULL) ? m_spatial->npars() : 0;
    int n_spectral = (m_spectral != NULL) ? m_spectral->npars() : 0;
    int n_temporal = (m_temporal != NULL) ? m_temporal->npars() : 0;
    m_npars        = n_spatial + n_spectral + n_temporal;

    // Continue only if there are parameters
    if (m_npars > 0) {

        // Allocate parameter pointers
        m_par = new GModelPar*[m_npars];

        // Initialise pointer on pointer array
        GModelPar** ptr = m_par;

        // Gather spatial parameter pointers
        for (int i = 0; i < n_spatial; ++i)
            *ptr++ = m_spatial->par(i);

        // Gather spectral parameters
        for (int i = 0; i < n_spectral; ++i)
            *ptr++ = m_spectral->par(i);

        // Gather temporal parameters
        for (int i = 0; i < n_temporal; ++i)
            *ptr++ = m_temporal->par(i);

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
        throw GException::model_invalid_spatial(G_XML_SPATIAL, type);

    // Type=ConstantValue
    else if (type == "ConstantValue")
        throw GException::model_invalid_spatial(G_XML_SPATIAL, type);

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
        throw GException::model_invalid_spectral(G_XML_SPECTRAL, type);

    // Unknown spectral model type
    else
        throw GException::model_invalid_spectral(G_XML_SPECTRAL, type);

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] obsDir Observed photon direction.
 * @param[in] obsEng Observed photon energy.
 * @param[in] obsTime Observed photon arrival time.
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @param[in] rsp Instrument response function.
 * @param[in] pnt Instrument pointing direction.
 * @parem[in] grad Compute also model gradients (default=false).
 *
 * This method evaluates the factorized source model at a given set of
 * parameters and multipies it with the corresponding value of the IRF.
 * Optionally, the parameter gradients are also evaluated and multiplied
 * with the corresponding IRF value. The result of this method needs to be
 * integrated over srcDir, srcEng, and srcTime to get the probability of
 * measuring an event.
 ***************************************************************************/
double GModel::fct(const GInstDir& obsDir, const GEnergy& obsEng,
                   const GTime& obsTime, const GSkyDir& srcDir,
                   const GEnergy& srcEng, const GTime& srcTime,
                   const GResponse& rsp, const GPointing& pnt, bool grad)
{
    // Initialise value
    double value = 0.0;

    // Compute IRF
    double irf = (&rsp)->irf(obsDir, obsEng, obsTime, srcDir, srcEng, srcTime, pnt);

    // Compute source model
    double source = 1.0;

    // Evaluate model in all components
    if (grad) {

        // Evaluate source model
        if (m_spatial  != NULL) source *= m_spatial->eval_gradients(srcDir);
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
        if (m_spatial  != NULL) source *= m_spatial->eval(srcDir);
        if (m_spectral != NULL) source *= m_spectral->eval(srcEng);
        if (m_temporal != NULL) source *= m_temporal->eval(srcTime);

        // Set value
        value = source * irf;
    }

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Perform integration over spatial component
 *
 * @param[in] obsDir Observed photon direction.
 * @param[in] obsEng Observed photon energy.
 * @param[in] obsTime Observed photon arrival time.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @param[in] rsp Instrument response function.
 * @param[in] pnt Instrument pointing direction.
 * @parem[in] grad Compute also model gradients (default=false).
 *
 * This method integrates the source model over the spatial component. If
 * the model does not explicitely depend on the sky direction (e.g. because
 * the model is a point source) no spatial integration is needed.
 *
 * @todo Needs implementation of spatial integration. Spatial integration
 * gets the integration region from the spatial model.
 ***************************************************************************/
double GModel::spatial(const GInstDir& obsDir, const GEnergy& obsEng,
                       const GTime& obsTime,
                       const GEnergy& srcEng, const GTime& srcTime,
                       const GResponse& rsp, const GPointing& pnt, bool grad)
{
    // Initialise result
    double value = 0.0;

    // Continue only if the gamma-ray source model has a spatial component
    if (m_spatial != NULL) {

        // Case A: Model is a point source
        if (m_spatial->isptsource()) {

            // Build sky direction from point source parameters
            GSkyDir srcDir;
            srcDir.radec_deg(((GModelSpatialPtsrc*)m_spatial)->ra(),
                             ((GModelSpatialPtsrc*)m_spatial)->dec());

            // Get function value at that position
            value = fct(obsDir, obsEng, obsTime, srcDir, srcEng, srcTime, rsp, pnt, grad);

        } // endif: Model was a point source

        // Case B: Model is not a point source
        else {

            // Dump warning that integration is not yet implemented
            std::cout << "WARNING: GModel::spatial:"
                      << " Sky integration not implemented." << std::endl;

        } // endelse: Model was not a point source

    } // endif: Gamma-ray source model had a spatial component

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Perform integration over spectral component
 *
 * @param[in] obsDir Observed photon direction.
 * @param[in] obsEng Observed photon energy.
 * @param[in] obsTime Observed photon arrival time.
 * @param[in] srcTime True photon arrival time.
 * @param[in] rsp Instrument response function.
 * @param[in] pnt Instrument pointing direction.
 * @parem[in] grad Compute also model gradients (default=false).
 *
 * This method integrates the source model over the spectral component. If
 * the response function has no energy dispersion then no spectral
 * integration is needed and the observed photon energy is identical to the
 * true photon energy.
 *
 * @todo Needs implementation of spectral integration.
 ***************************************************************************/
double GModel::spectral(const GInstDir& obsDir, const GEnergy& obsEng,
                        const GTime& obsTime, const GTime& srcTime,
                        const GResponse& rsp, const GPointing& pnt, bool grad)
{
    // Initialise result
    double value = 0.0;

    // Determine if energy integration is needed
    bool integrate = rsp.hasedisp();

    // Case A: Integraion
    if (integrate) {
        std::cout << "GModel::spectral: Integration not implemented." << std::endl;
    }

    // Case B: No integration (assume no energy dispersion)
    else
        value = spatial(obsDir, obsEng, obsTime, obsEng, srcTime, rsp, pnt, grad);

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Perform integration over temporal component
 *
 * @param[in] obsDir Observed photon direction.
 * @param[in] obsEng Observed photon energy.
 * @param[in] obsTime Observed photon arrival time.
 * @param[in] rsp Instrument response function.
 * @param[in] pnt Instrument pointing direction.
 * @parem[in] grad Compute also model gradients (default=false).
 *
 * This method integrates the source model over the temporal component. If
 * the response function has no time dispersion then no temporal integration
 * is needed and the observed photon arrival time is identical to the true
 * photon arrival time.
 *
 * @todo Needs implementation of temporal integration.
 ***************************************************************************/
double GModel::temporal(const GInstDir& obsDir, const GEnergy& obsEng,
                        const GTime& obsTime,
                        const GResponse& rsp, const GPointing& pnt, bool grad)
{
    // Initialise result
    double value = 0.0;

    // Determine if time integration is needed
    bool integrate = rsp.hastdisp();

    // Case A: Integraion
    if (integrate) {
        std::cout << "GModel::temporal: Integration not implemented." << std::endl;
    }

    // Case B: No integration (assume no time dispersion)
    else
        value = spectral(obsDir, obsEng, obsTime, obsTime, rsp, pnt, grad);

    // Return value
    return value;
}


/*==========================================================================
 =                                                                         =
 =                               GModel friends                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put model in output stream
 *
 * @param[in] os Output stream into which the model will be dumped
 * @param[in] model Model to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GModel& model)
{
    // Determine number of parameters per type
   int n_spatial  = model.m_spatial->npars();
   int n_spectral = model.m_spectral->npars();
   int n_temporal = model.m_temporal->npars();

    // Put model in stream
    os << "=== GModel ===" << std::endl;
    os << " Name ......................: " << model.m_name << std::endl;
    os << " Number of parameters ......: " << model.m_npars << std::endl;
    os << " Number of spatial par's ...: " << n_spatial;
    int i;
    for (i = 0; i < n_spatial; ++i) {
        os << std::endl;
        os << " Parameter .................: " << *(model.m_par[i]);
    }
    os << std::endl << " Number of spectral par's ..: " << n_spatial;
    for (; i < n_spatial+n_spectral; ++i) {
        os << std::endl;
        os << " Parameter .................: " << *(model.m_par[i]);
    }
    os << std::endl << " Number of temporal par's ..: " << n_temporal;
    for (; i < n_spatial+n_spectral+n_temporal; ++i) {
        os << std::endl;
        os << " Parameter .................: " << *(model.m_par[i]);
    }    

    // Return output stream
    return os;
}
