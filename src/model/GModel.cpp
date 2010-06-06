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

/* __ Method name definitions ____________________________________________ */
#define G_PAR                                        "GModel::par(int) const"
#define G_SET_POINTERS                               "GModel::set_pointers()"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                       GModel constructors/destructors                   =
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
 =                             GModel operators                            =
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
 =                          GModel public methods                          =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Get pointer to model parameter
 *
 * @param[in] index Parameter index.
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
 * @brief Return source model value
 *
 * @param[in] obsDir Observed photon direction (to be removed).
 * @param[in] obsEng Observed photon energy (to be removed).
 * @param[in] obsTime Observed photon arrival time (to be removed).
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @param[in] rsp Instrument response function (to be removed).
 * @param[in] pnt Instrument pointing direction (to be removed).
 *
 * This method evaluates the factorized source model at a given set of
 * parameters.
 *
 * @todo obsDir, obsEng, obsTime, rsp and pnt arguments are to be removed.
 ***************************************************************************/
double GModel::value(const GInstDir& obsDir, const GEnergy& obsEng,
                     const GTime& obsTime, const GSkyDir& srcDir,
                     const GEnergy& srcEng, const GTime& srcTime,
                     const GResponse& rsp, const GPointing& pnt)
{
    // Initialise source model
    double source = 1.0;

    // Evaluate source model
    if (m_spatial  != NULL)
        source *= m_spatial->eval(obsDir, srcDir, srcEng, srcTime, rsp, pnt);
    if (m_spectral != NULL)
        source *= m_spectral->eval(obsEng, srcDir, srcEng, srcTime, rsp, pnt);
    if (m_temporal != NULL)
        source *= m_temporal->eval(obsTime, srcDir, srcEng, srcTime, rsp, pnt);

    // Return
    return source;
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
        if (m_par == NULL)
            throw GException::mem_alloc(G_SET_POINTERS, m_npars);
    
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
        if (m_spatial  != NULL)
            source *= m_spatial->eval_gradients(obsDir, srcDir, srcEng, srcTime, rsp, pnt);
        if (m_spectral != NULL)
            source *= m_spectral->eval_gradients(obsEng, srcDir, srcEng, srcTime, rsp, pnt);
        if (m_temporal != NULL)
            source *= m_temporal->eval_gradients(obsTime, srcDir, srcEng, srcTime, rsp, pnt);

        // Set value
        value = source * irf;
    
        // Set gradients
        for (int i = 0; i < m_npars; ++i)
            m_par[i]->gradient(m_par[i]->gradient() * irf);
    
    }
    else {
    
        // Evaluate source model
        if (m_spatial  != NULL)
            source *= m_spatial->eval(obsDir, srcDir, srcEng, srcTime, rsp, pnt);
        if (m_spectral != NULL)
            source *= m_spectral->eval(obsEng, srcDir, srcEng, srcTime, rsp, pnt);
        if (m_temporal != NULL)
            source *= m_temporal->eval(obsTime, srcDir, srcEng, srcTime, rsp, pnt);

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
 * @todo Quick and dirty source location extraction from GModelSpatialPtsrc.
 ***************************************************************************/
double GModel::spatial(const GInstDir& obsDir, const GEnergy& obsEng,
                       const GTime& obsTime,
                       const GEnergy& srcEng, const GTime& srcTime,
                       const GResponse& rsp, const GPointing& pnt, bool grad)
{
    // Initialise result
    double value = 0.0;

    // Determine if integration is needed
    bool integrate  = (m_spatial != NULL) ? m_spatial->depdir() : false;

    // Case A: Integraion
    if (integrate) {
        std::cout << "GModel::spatial: Integration not implemented." << std::endl;
    }
    
    // Case B: No integration, then extract point source position from model
    else {        
        GSkyDir srcDir;
        srcDir.radec_deg(((GModelSpatialPtsrc*)m_spatial)->ra(),
                         ((GModelSpatialPtsrc*)m_spatial)->dec());
        value = fct(obsDir, obsEng, obsTime, srcDir, srcEng, srcTime, rsp, pnt, grad);
    }
    
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


/*==========================================================================
 =                                                                         =
 =                      Other functions used by GModel                     =
 =                                                                         =
 ==========================================================================*/
