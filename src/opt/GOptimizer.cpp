/***************************************************************************
 *           GOptimizer.cpp  -  Abstract base class for optimizer          *
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
 * @file GOptimizer.cpp
 * @brief GOptimizer abstract base class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GOptimizer.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_CONSTRUCTOR        "GOptimizer::GOptimizer(const GModels&)"
#define G_PAR                                "GOptimizer::par(int) const"
#define G_COPY_MEMBERS          "GOptimizer::copy_members(const GOptimizer&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                    GOptimizer constructors/destructors                  =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GOptimizer::GOptimizer(void)
{
    // Initialise private members for clean destruction
    init_members();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Construct object from function and parameters
 *
 * @param[in] fct Optimizer function.
 * @param[in] pars Optimizer parameters.
 ***************************************************************************/
GOptimizer::GOptimizer(const GOptimizerFunction& fct, const GOptimizerPars &pars)
{
    // Initialise private members for clean destruction
    init_members();
    
    // Store optimizer function and parameters
    m_fct  = &fct;
    m_pars = pars;
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Construct object from function and models
 *
 * @param[in] fct Optimizer function.
 * @param[in] models Optimizer parameters.
 ***************************************************************************/
GOptimizer::GOptimizer(const GOptimizerFunction& fct, const GModels &models)
{
    // Initialise private members for clean destruction
    init_members();

    // Store optimizer function and parameters
    m_fct  = &fct;
    m_pars = GOptimizerPars(models);
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] opt Optimizer from which the instance should be built.
 ***************************************************************************/
GOptimizer::GOptimizer(const GOptimizer& opt)
{ 
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(opt);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GOptimizer::~GOptimizer()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                           GOptimizer operators                          =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] opt Optimizer to be assigned.
 ***************************************************************************/
GOptimizer& GOptimizer::operator= (const GOptimizer& opt)
{ 
    // Execute only if object is not identical
    if (this != &opt) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(opt);

    } // endif: object was not identical
  
    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                         GOptimizer public methods                       =
 =                                                                         =
 ==========================================================================*/

/*==========================================================================
 =                                                                         =
 =                         GOptimizer private methods                      =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GOptimizer::init_members(void)
{
    // Initialise members
    m_fct  = NULL;
    m_pars = GOptimizerPars();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] opt GOptimizer members to be copied.
 ***************************************************************************/
void GOptimizer::copy_members(const GOptimizer& opt)
{
    // Copy attributes
    m_fct  = opt.m_fct;
    m_pars = opt.m_pars;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GOptimizer::free_members(void)
{
    // Free memory

    // Signal free pointers
  
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             GOptimizer friends                          =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put optimizer in output stream
 *
 * @param[in] os Output stream into which the optimizer will be dumped
 * @param[in] pars Parameters to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GOptimizer& opt)
{
    // Put optimizer in stream
    os << "=== GOptimizer ===" << std::endl;

    // Return output stream
    return os;
}


/*==========================================================================
 =                                                                         =
 =                    Other functions used by GOptimizer                   =
 =                                                                         =
 ==========================================================================*/
