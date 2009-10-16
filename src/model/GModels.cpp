/***************************************************************************
 *                  GModels.cpp  -  Model container class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
/**
 * @file GModels.cpp
 * @brief GModels class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GModels.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ADD                             "GModels::add(const GModel& model)"
#define G_PAR                                       "GModels::par(int) const"
#define G_COPY_MEMBERS                "GModels::copy_members(const GModels&)"
#define G_SET_POINTERS                              "GModels::set_pointers()"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                      GModels constructors/destructors                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
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
 * @param[in] models Models from which the instance should be built.
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
 * @brief Destructor
 ***************************************************************************/
GModels::~GModels()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                            GModels operators                            =
 =                                                                         =
 ==========================================================================*/

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
 =                         GModels public methods                          =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Add model to container
 *
 * @param[in] model Model to be added.
 ***************************************************************************/
void GModels::add(const GModel& model)
{
	// Allocate fresh memory for models
    GModel* new_model = new GModel[m_elements+1];
	if (new_model == NULL)
		throw GException::mem_alloc(G_ADD, m_elements+1);

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


/*==========================================================================
 =                                                                         =
 =                         GModels private methods                         =
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
        if (m_model == NULL)
            throw GException::mem_alloc(G_COPY_MEMBERS, m_elements);
        
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
 ***************************************************************************/
void GModels::set_pointers(void)
{
    // Delete old pointers (if they exist)
    if (m_par != NULL) delete [] m_par;
    m_par = NULL;
    
    // Determine total number of parameters
    m_npars = 0;
    for (int k = 0; k < m_elements; ++k)
        m_npars += m_model[k].npars();
    
    // Continue only if there are parameters
    if (m_npars > 0) {
    
        // Allocate parameter pointers
        m_par = new GModelPar*[m_npars];
        if (m_par == NULL)
            throw GException::mem_alloc(G_SET_POINTERS, m_npars);
    
        // Initialise pointer on pointer array
        GModelPar** ptr = m_par;
        
        // Loop over all models
        for (int k = 0; k < m_elements; ++k) {

            // Determine the number of parameters of each type
            int n_spatial  = (m_model[k].m_spatial  != NULL) 
                             ? m_model[k].m_spatial->npars() : 0;
            int n_spectral = (m_model[k].m_spectral != NULL) 
                             ? m_model[k].m_spectral->npars() : 0;
            int n_temporal = (m_model[k].m_temporal != NULL) 
                             ? m_model[k].m_temporal->npars() : 0;
        
            // Gather spatial parameter pointers
            for (int i = 0; i < n_spatial; ++i)
                *ptr++ = m_model[k].m_spatial->par(i);
        
            // Gather spectral parameters
            for (int i = 0; i < n_spectral; ++i)
                *ptr++ = m_model[k].m_spectral->par(i);

            // Gather temporal parameters
            for (int i = 0; i < n_temporal; ++i)
                *ptr++ = m_model[k].m_temporal->par(i);
        
        } // endfor: looped over all models
        
    } // endif: there were model parameters
    
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                              GModels friends                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put models in output stream
 *
 * @param[in] os Output stream into which the model will be dumped
 * @param[in] model Model to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GModels& models)
{
    // Allocate filler
    std::string filler = " ..............";
    
    // Put model in stream
    os << "=== GModels ===" << std::endl;
    os << " Number of models ..........: " << models.m_elements << std::endl;
    os << " Number of parameters ......: " << models.m_npars;
    int i = 0;
    for (int k = 0; k < models.m_elements; ++k) {
        os << std::endl << " Model name ................: " << models.m_model[k].name();
        for (int j = 0; j < models.m_model[k].npars(); ++j, ++i) {
            os << std::endl;
            if (i == 10) filler = " .............";
            if (i == 100) filler = " ............";
            if (i == 1000) filler = " ...........";
            os << "  Parameter " << i << filler << ": " << *(models.m_par[i]);
        }
    }

    // Return output stream
    return os;
}


/*==========================================================================
 =                                                                         =
 =                      Other functions used by GModels                     =
 =                                                                         =
 ==========================================================================*/
