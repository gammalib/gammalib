/***************************************************************************
 *            GOptimizerLM.cpp  -  Levenberg Marquardt optimizer           *
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
 * @file GOptimizerLM.cpp
 * @brief GOptimizerLM base class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GOptimizerLM.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
#define G_DEBUG_OPT 1 //!< Perform optimizer debugging (0=no, 1=yes)


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
 * @param[in] p Parameters to be optimised.
 ***************************************************************************/
GOptimizerPars& GOptimizerLM::operator() (GOptimizerFunction& fct, GOptimizerPars& p)
{
    // Initalise output parameters with input parameters
    GOptimizerPars* pars = new GOptimizerPars(p);
    
    // Perform LM optimization
    optimize(&fct, pars);
    
    // Return
    return *pars;
}


/***********************************************************************//**
 * @brief Optimization operator
 *
 * @param[in] fct Optimization function.
 * @param[in] m Model parameters to be optimised.
 ***************************************************************************/
GModels& GOptimizerLM::operator() (GOptimizerFunction& fct, GModels& m)
{
    // Initalise output parameters with input parameters
    GModels* models = new GModels(m);
    
    // Perform LM optimization
    optimize(&fct, models);

    // Return
    return *models;
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
    // Initialise optimizer parameters
    m_lambda_start = 1.0e-3;
    m_lambda_inc   = 10000.0;
    m_lambda_dec   = 0.1;
    m_eps          = 1.0e-6;
    m_max_iter     = 100;
    m_max_stall    = 10;

    // Initialise optimizer values
    m_lambda = m_lambda_start;
    m_value  = 0.0;
    m_status = 0;
    m_iter   = 0;

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


/***********************************************************************//**
 * @brief Perform LM optimization
 *
 * @param[in] fct Poiner to optimization function.
 * @param[in] par Pointer to parameters to be optimised.
 *
 * @todo Implemenet logger to allow for optimizer logging in applications.
 * So far only use std::cout (definition enabled).
 ***************************************************************************/
void GOptimizerLM::optimize(GOptimizerFunction* fct, GOptimizerPars* pars)
{
    // Single loop for common exit point
    do {
    
        // Initialise optimization parameters
        m_lambda = m_lambda_start;
        m_status = 0;
        
        // Initial evaluation
        fct->eval(*pars);
        
        // Save parameters
        m_value = *(fct->value());
        
        // Save initial statistics and lambda values
        double value_old  = m_value;
        double lambda_old = m_lambda;
        int    lambda_inc = 0;
        
        // Dump
        #if G_DEBUG_OPT
        std::cout << "Initial iteration: func=" << m_value << ", Lambda=" 
                  << m_lambda << std::endl;
        #endif
        
        // Iterative fitting
        for (m_iter = 0; m_iter < m_max_iter; ++m_iter) {
        
            // Perform one iteration
            iteration(fct, pars);
            
            // Compute function improvement (>0 means decrease)
            double delta = value_old - m_value;
            
            // Dump
            #if G_DEBUG_OPT
            std::cout << "Iteration " << m_iter+1 << ": func=" 
                      << m_value << ", Lambda=" << m_lambda
                      << ", delta=" << delta << std::endl;
            #endif
            
            // Reset lambda increment if we had success
            if (m_lambda < lambda_old)
                lambda_inc = 0;

            // If function increased while lambda did not increase then stop
            // iterations
            if ((m_lambda <= lambda_old) && (delta < 0.0))
                break;

            // Stop if convergence was reached
            if ((m_lambda <= lambda_old) && (delta < m_eps))
                break;

            // Monitor the number of subsequent increases of lambda and stop if
            // the number of increases exceeds threshold
            lambda_inc = (m_lambda > lambda_old) ? lambda_inc + 1 : 0;
            if (lambda_inc > m_max_stall) {
                m_status = G_LM_STALLED;
                break;
            }

            // Bookkeeping of actual result (we always store the last lambda to
            // detect turn arounds in the lambda tendency; however we always keep
            // the best function value)
            lambda_old = m_lambda;
            if (delta > 0.0)
                value_old = m_value;

        } // endfor: iterations
    
    } while (0); // endwhile: main loop
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Perform one LM iteration
 *
 * @param[in] fct Poiner to optimization function.
 * @param[in] par Pointer to parameters to be optimised.
 ***************************************************************************/
void GOptimizerLM::iteration(GOptimizerFunction* fct, GOptimizerPars* pars)
{
    // Single loop for common exit point
    do {
        // Initialise iteration parameters
        int            npars = pars->npars();
        GVector*       grad  = fct->gradient();
        GSparseMatrix* covar = fct->covar();
    
        // Save function value, gradient and covariance matrix
        double         save_value = m_value;
        GVector        save_grad  = GVector(*grad);
        GSparseMatrix  save_covar = GSparseMatrix(*covar);
        
        // Save parameter values in vector
        GVector save_pars(npars);
        for (int ipar = 0; ipar < npars; ++ipar)
            save_pars(ipar) = pars->par(ipar)->value();
        
        // Setup matrix and vector for covariance computation
        for (int ipar = 0; ipar < npars; ++ipar) {
            (*covar)(ipar,ipar) *= (1.0 + m_lambda);
            (*grad)(ipar)        = -(*grad)(ipar);
        }
        
        // Solve: covar * X = grad. Handle matrix problems
        try {
            (*covar).cholesky_decompose(1);
            *grad = (*covar).cholesky_solver(*grad);
        }
        catch (GException::matrix_zero &e) {
            m_status = G_LM_SINGULAR;
            std::cout << "GOptimizerLM::iteration: all curvature matrix"
                      << " elements are 0." << std::endl;
            continue;
        }
        catch (GException::matrix_not_pos_definite &e) {
            m_status = G_LM_NOT_POSTIVE_DEFINITE;
            std::cout << "GOptimizerLM::iteration: Curvature matrix not"
                      << " positive definite." << std::endl;
            continue;
        }
        catch (std::exception &e) {
            throw;
        }
        
        // Derive new parameter vector
        double step = 1.0; // NOTE: STEP SIZE ADJUSTMENT TBW

        // Derive new parameter vector
        for (int ipar = 0; ipar < npars; ++ipar) {
        
            // Get actual parameter value, limits and scale
            double p     = pars->par(ipar)->value();
            double p_min = pars->par(ipar)->min();
            double p_max = pars->par(ipar)->max();
            double scale = pars->par(ipar)->scale();
        
            // Compute new parameter value
            p += (*grad)(ipar) * scale * step;
            
            // Constrain parameter to within the valid range
            if (pars->par(ipar)->hasmin() && p < p_min) {
                std::cout << "Parameter " << ipar << " hits minimum: " << p << " < "
                          << p_min << std::endl;
                p = p_min;
            }
            if (pars->par(ipar)->hasmax() && p > p_max) {
                std::cout << "Parameter " << ipar << " hits maximum: " << p << " > "
                          << p_max << std::endl;
                p = p_max;
            }
            
            // Set parameter value
            pars->par(ipar)->value(p);
            
        } // endfor: computed new parameter vector
        
        // Evaluate function at new parameters
        fct->eval(*pars);
        
        // Fetch new pointers since eval will allocate new memory
        grad  = fct->gradient();
        covar = fct->covar();
        
        // Retrieve new function value
        m_value = *(fct->value());

        // Determine how many parameters have changed
        int par_change = 0;
        for (int ipar = 0; ipar < npars; ++ipar) {
            if (pars->par(ipar)->value() != save_pars(ipar))
                par_change++;
        }

        // If the function has decreased then accept the new solution and decrease
        // lambda ...
        double delta = save_value - m_value;
        if (delta > 0.0)
            m_lambda *= m_lambda_dec;
            
        // ... if function is identical then accept new solution. If the parameters 
        // have changed then increase lambda, otherwise decrease lambda
        else if (delta == 0.0) {
            if (par_change)
                m_lambda *= m_lambda_inc;
            else
                m_lambda *= m_lambda_dec;
        }
        
        // ... if function worsened slightly then accept new solution and increase
        // lambda
        else if (delta > 1.0e-6)
            m_lambda *= m_lambda_inc;
        
        // ... otherwise,  if the statistics did not improve then use old parameters
        // and increase lamdba. Restore also the best statistics value that was 
        // reached so far, the gradient vector and the curve matrix.
        else {
            m_lambda *= m_lambda_inc;
            m_value   = save_value;
            *grad     = save_grad;
            *covar    = save_covar;
            for (int ipar = 0; ipar < npars; ++ipar)
                pars->par(ipar)->value(save_pars(ipar));
        }

    } while (0); // endwhile: main loop
    
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
