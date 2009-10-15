/***************************************************************************
 *            GOptimizerLM.cpp  -  Levenberg Marquardt optimizer           *
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
 * @file GOptimizerLM.cpp
 * @brief GOptimizerLM base class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GOptimizerLM.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                   GOptimizerLM constructors/destructors                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GOptimizerLM::GOptimizerLM(void) : GOptimizer()
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
/*
GOptimizerLM::GOptimizerLM(const GOptimizerFunction& fct, const GOptimizerPars &pars) :
  GOptimizer(fct, pars)
{
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}
*/

/***********************************************************************//**
 * @brief Construct object from function and models
 *
 * @param[in] fct Optimizer function.
 * @param[in] models Optimizer parameters.
 ***************************************************************************/
/*
GOptimizerLM::GOptimizerLM(const GOptimizerFunction& fct, const GModels &models) :
  GOptimizer(fct, models)
{
    // Initialise private members for clean destruction
    init_members();
  
    // Return
    return;
}
*/

/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] opt Optimizer from which the instance should be built.
 ***************************************************************************/
GOptimizerLM::GOptimizerLM(const GOptimizerLM& opt) : GOptimizer(opt)
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
GOptimizerLM::~GOptimizerLM()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                          GOptimizerLM operators                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] opt Optimizer to be assigned.
 ***************************************************************************/
GOptimizerLM& GOptimizerLM::operator= (const GOptimizerLM& opt)
{ 
    // Execute only if object is not identical
    if (this != &opt) {

        // Copy base class members
        this->GOptimizer::operator=(opt);

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


/***********************************************************************//**
 * @brief Optimization operator
 *
 * @param[in] fct Optimization function.
 * @param[in] pars Parameters to be optimised.
 ***************************************************************************/
GOptimizerPars& GOptimizerLM::operator() (GOptimizerFunction& fct, GOptimizerPars& p)
{
    // Initalise output parameters with input parameters
    GOptimizerPars* pars = new GOptimizerPars(p);
    
    // DUMMY
    pars->par(0)->value(10.0);
    pars->par(0)->error(1.0);

    // Return
    return *pars;
}


/*==========================================================================
 =                                                                         =
 =                        GOptimizerLM public methods                      =
 =                                                                         =
 ==========================================================================*/

/*==========================================================================
 =                                                                         =
 =                        GOptimizerLM private methods                     =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GOptimizerLM::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] opt GOptimizerLM members to be copied.
 ***************************************************************************/
void GOptimizerLM::copy_members(const GOptimizerLM& opt)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GOptimizerLM::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                            GOptimizerLM friends                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put optimizer in output stream
 *
 * @param[in] os Output stream into which the optimizer will be dumped
 * @param[in] pars Parameters to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GOptimizerLM& opt)
{
    // Put optimizer in stream
    os << "=== GOptimizerLM ===" << std::endl;

    // Return output stream
    return os;
}


/*==========================================================================
 =                                                                         =
 =                    Other functions used by GOptimizer                   =
 =                                                                         =
 ==========================================================================*/
