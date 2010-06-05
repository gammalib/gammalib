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
 * @brief Evaluate function
 *
 * @param[in] obsDir Observed photon direction.
 * @param[in] obsEng Observed photon energy.
 * @param[in] obsTime Observed photon arrival time.
 * @param[in] rsp Instrument response function.
 * @param[in] pnt Instrument pointing direction.
 *
 * This method implements a source model evaluation assuming a factorisation
 * of the source term into a spatial, a spectral and a temporal component.
 *
 * @todo General method needs to be implemented.
 ***************************************************************************/
double GModel::eval(const GInstDir& obsDir, const GEnergy& obsEng,
                    const GTime& obsTime, const GResponse& rsp,
                    const GPointing& pnt)
{
    // Initialise value
    double value = 1.0;

    // Integral over source direction, energy and time
    //for (srcTime = ...
    //for (srcEng = ...
    //for (srcDir = ...
    GTime   srcTime = obsTime;    // Assume no time dispersion
    GEnergy srcEng  = obsEng;     // Assume no energy dispersion
    GSkyDir srcDir;               // Assume nothing
    
    // Evaluate model in all components
    if (m_spatial  != NULL)
        value *= m_spatial->eval(obsDir, srcDir, srcEng, srcTime, rsp, pnt);
    if (m_spectral != NULL)
        value *= m_spectral->eval(obsEng, srcDir, srcEng, srcTime, rsp, pnt);
    if (m_temporal != NULL)
        value *= m_temporal->eval(obsTime);
    
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
 * This method implements a source model evaluation assuming a factorisation
 * of the source term into a spatial, a spectral and a temporal component.
 *
 * @todo General method needs to be implemented.
 ***************************************************************************/
double GModel::eval_gradients(const GInstDir& obsDir, const GEnergy& obsEng,
                              const GTime& obsTime, const GResponse& rsp,
                              const GPointing& pnt)
{
    // Initialise value
    double value = 1.0;

    // Integral over source direction, energy and time
    //for (srcTime = ...
    //for (srcEng = ...
    //for (srcDir = ...
    GTime   srcTime = obsTime;    // Assume no time dispersion
    GEnergy srcEng  = obsEng;     // Assume no energy dispersion
    GSkyDir srcDir;               // Assume nothing
    
    // Evaluate model and gradients in all components
    if (m_spatial  != NULL) 
        value *= m_spatial->eval_gradients(obsDir, srcDir, srcEng, srcTime, rsp, pnt);
    if (m_spectral != NULL)
        value *= m_spectral->eval_gradients(obsEng, srcDir, srcEng, srcTime, rsp, pnt);
    if (m_temporal != NULL)
        value *= m_temporal->eval_gradients(obsTime);
    
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
