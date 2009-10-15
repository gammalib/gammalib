/***************************************************************************
 *          GData_optimizer.cpp  -  Optimizer class of data class          *
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
 * @file GData_optimizer.cpp
 * @brief GData::optimizer class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GData.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Prototypes _________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                GData::optimizer constructors/destructors                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
***************************************************************************/
GData::optimizer::optimizer(void) : GOptimizerFunction()
{
    // Initialise iterator
    init_members();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor based on specific GData object
 ***************************************************************************/
GData::optimizer::optimizer(GData *data) : GOptimizerFunction()
{
    // Initialise iterator
    init_members();

    // Set data object
    m_data = data;
    
    // Return
    return;
}



/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param data GData instance which should be used for construction
 ***************************************************************************/
GData::optimizer::optimizer(const optimizer& fct) : GOptimizerFunction(fct)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(fct);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GData::optimizer::~optimizer()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                       GData::optimizer operators                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] data GData instance to be assigned
 ***************************************************************************/
GData::optimizer& GData::optimizer::operator= (const optimizer& fct)
{
    // Execute only if object is not identical
    if (this != &fct) {

        // Copy base class members
        this->GOptimizerFunction::operator=(fct);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(fct);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                      GData::optimizer public methods                    =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Evaluate function
 ***************************************************************************/
void GData::optimizer::eval(const GOptimizerPars& pars) 
{
    // DUMMY
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Return function value
 ***************************************************************************/
double* GData::optimizer::value(void) 
{
    // DUMMY
    
    // Return
    return &m_value;
}


/***********************************************************************//**
 * @brief Return function gradient
 ***************************************************************************/
GVector* GData::optimizer::gradient(void) 
{
    // DUMMY
    
    // Return
    return m_vector;
}


/***********************************************************************//**
 * @brief Return function covariance matrix
 ***************************************************************************/
GSparseMatrix* GData::optimizer::covar(void) 
{
    // DUMMY
    
    // Return
    return m_covar;
}


/*==========================================================================
 =                                                                         =
 =                     GData::optimizer private methods                    =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GData::optimizer::init_members(void)
{
    // Initialise members
    m_value  = 0.0;
    m_vector = NULL;
    m_covar  = NULL;
    m_data   = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] fct GData::optimizer members to be copied.
 ***************************************************************************/
void GData::optimizer::copy_members(const optimizer& fct)
{
    // Copy attributes
    m_value = fct.m_value;
    
    // Copy gradient if it exists
    if (fct.m_vector != NULL)
        m_vector = new GVector(*fct.m_vector);

    // Copy covariance matrix if it exists
    if (fct.m_covar != NULL)
        m_covar = new GSparseMatrix(*fct.m_covar);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GData::optimizer::free_members(void)
{
    // Free members
	if (m_vector != NULL) delete m_vector;
	if (m_covar  != NULL) delete m_covar;
    
    // Signal free pointers
    m_vector = NULL;
    m_covar  = NULL;

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                         GData::optimizer friends                        =
 =                                                                         =
 ==========================================================================*/

/*==========================================================================
 =                                                                         =
 =                  Other functions used by GData::optimizer               =
 =                                                                         =
 ==========================================================================*/
