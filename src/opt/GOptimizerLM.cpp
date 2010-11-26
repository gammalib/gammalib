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
#define G_DEBUG_OPT 0          //!< Perform optimizer debugging (0=no, 1=yes)


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GOptimizerLM::GOptimizerLM(void) : GOptimizer()
{
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor with logger
 *
 * @param[in] log Logger to use in optimizer.
 ***************************************************************************/
GOptimizerLM::GOptimizerLM(GLog& log) : GOptimizer()
{
    // Initialise private members for clean destruction
    init_members();

    // Set pointer to logger
    m_logger = &log;

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
 =                               Operators                                 =
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
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
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
    m_step_adjust  = true;

    // Initialise optimizer values
    m_lambda = m_lambda_start;
    m_value  = 0.0;
    m_status = 0;
    m_iter   = 0;

    // Initialise pointer to logger
    m_logger = NULL;

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
    // Copy attributes
    m_lambda_start = opt.m_lambda_start;
    m_lambda_inc   = opt.m_lambda_inc;
    m_lambda_dec   = opt.m_lambda_dec;
    m_eps          = opt.m_eps;
    m_max_iter     = opt.m_max_iter;
    m_max_stall    = opt.m_max_stall;
    m_step_adjust  = opt.m_step_adjust;
    m_lambda       = opt.m_lambda;
    m_value        = opt.m_value;
    m_status       = opt.m_status;
    m_iter         = opt.m_iter;
    m_logger       = opt.m_logger;

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
 * @param[in] fct Optimization function.
 * @param[in] pars Function parameters.
 ***************************************************************************/
void GOptimizerLM::optimize(GOptimizerFunction* fct, GOptimizerPars* pars)
{
    // Initialise optimization parameters
    int npars = pars->npars();
    m_lambda  = m_lambda_start;
    m_status  = 0;

    // Allocate temporary memory
    m_hit_boundary = new bool[npars];
    for (int i = 0; i < npars; ++i)
        m_hit_boundary[i] = false;
    
    // Initial evaluation
    fct->eval(*pars);

    // Save parameters
    m_value = *(fct->value());

    // Save initial statistics and lambda values
    double value_old  = m_value;
    double lambda_old = m_lambda;
    int    lambda_inc = 0;

    // Optionally write initial iteration into logger
    if (m_logger != NULL) {
        *m_logger << "Initial iteration: ";
        *m_logger << "func=" << m_value << ", ";
        *m_logger << "Lambda=" << m_lambda << std::endl;
    }
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

        // Optionally write iteration results into logger
        if (m_logger != NULL) {
            *m_logger << "Iteration " << m_iter+1 << ": ";
            *m_logger << "func=" << m_value << ", ";
            *m_logger << "Lambda=" << m_lambda << ", ";
            *m_logger << "delta=" << delta << std::endl;
        }
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

    // Compute parameter uncertainties
    errors(fct, pars);

    // Free working memory
    if (m_hit_boundary != NULL) delete [] m_hit_boundary;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Perform one LM iteration
 *
 * @param[in] fct Optimizer function.
 * @param[in] pars Function parameters.
 *
 * This method performs one LM iteration. Note that the method only acts on
 * the parameter value, i.e. it does not worry about the true scale of the
 * parameter. It calls the eval() method of the optimizer function which
 * is assumed to return the gradient and curvature matrix with respect to the
 * parameter values (and not the scaled true values).
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
            covar->cholesky_decompose(1);
            *grad = covar->cholesky_solver(*grad);
        }
        catch (GException::matrix_zero &e) {
            m_status = G_LM_SINGULAR;
            if (m_logger != NULL) {
                *m_logger << "GOptimizerLM::iteration: "
                          << "All curvature matrix elements are zero."
                          << std::endl;
            }
            continue;
        }
        catch (GException::matrix_not_pos_definite &e) {
            m_status = G_LM_NOT_POSTIVE_DEFINITE;
            if (m_logger != NULL) {
                *m_logger << "GOptimizerLM::iteration: "
                          << "Curvature matrix not positive definite."
                          << std::endl;
            }
            continue;
        }
        catch (std::exception &e) {
            throw;
        }

        // Get LM step size
        double step = step_size(grad, pars);
     
        // Derive new parameter vector
        for (int ipar = 0; ipar < npars; ++ipar) {

            // Get actual parameter value and limits
            double p     = pars->par(ipar)->value();
            double p_min = pars->par(ipar)->min();
            double p_max = pars->par(ipar)->max();

            // Compute new parameter value
            p += (*grad)(ipar) * step;

            // Constrain parameter to within the valid range
            if (pars->par(ipar)->hasmin() && p < p_min) {
                if (m_logger != NULL) {
                    *m_logger << "... parameter \"" << pars->par(ipar)->name();
                    *m_logger << "\" hits minimum: ";
                    *m_logger << p << " < " << p_min << std::endl;
                }
                p = p_min;
            }
            if (pars->par(ipar)->hasmax() && p > p_max) {
                if (m_logger != NULL) {
                    *m_logger << "... parameter \"" << pars->par(ipar)->name();
                    *m_logger << "\" hits maximum: ";
                    *m_logger << p << " > " << p_max << std::endl;
                }
                p = p_max;
            }

            // Set new parameter value
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


/***********************************************************************//**
 * @brief Return LM step size
 *
 * @param[in] grad Function gradient.
 * @param[in] pars Function parameters.
 *
 * Determine the size of the LM step. By default a step size of 1 is taken.
 * If m_step_adjust=true then the step size will be estimated so that the
 * next parameter vector should stay within the parameter boundaries
 * (provided that boundaries exist).
 ***************************************************************************/
double GOptimizerLM::step_size(GVector* grad, GOptimizerPars* pars)
{
    // Initialise step size
    double step = 1.0;

    // Check if we should reduce the step size
    if (m_step_adjust) {

        // Initialise the parameter index that constrains most the fit
        int ipar_bnd = -1;

        // Loop over all parameters
        for (int ipar = 0; ipar < pars->npars(); ++ipar) {

            // Get parameter attributes
            double p     = pars->par(ipar)->value();
            double p_min = pars->par(ipar)->min();
            double p_max = pars->par(ipar)->max();
            double delta = (*grad)(ipar);

            // Check if a parameter minimum requires a reduced step size
            if (pars->par(ipar)->hasmin()) {
                double step_min = (delta < 0.0) ? (p_min - p)/delta : 1.0;
                if (step_min > 0.0) {
                    if (step_min < step) {
                        if (!m_hit_boundary[ipar]) {
                            ipar_bnd = ipar;
                            step     = step_min;
                        }
                    }
                    else if (m_hit_boundary[ipar]) {
                        m_hit_boundary[ipar] = false;
                        if (m_logger != NULL) {
                            *m_logger << "... parameter \"";
                            *m_logger << pars->par(ipar)->name();
                            *m_logger << "\" does not drive optimization step anymore.";
                            *m_logger << std::endl;
                        }
                    }
                }
            }

            // Check if a parameter maximum requires a reduced step size
            if (pars->par(ipar)->hasmax()) {
                double step_max = (delta > 0.0) ? (p_max - p)/delta : 1.0;
                if (step_max > 0.0) {
                    if (step_max < step) {
                        if (!m_hit_boundary[ipar]) {
                            ipar_bnd = ipar;
                            step     = step_max;
                        }
                    }
                    else if (m_hit_boundary[ipar]) {
                        m_hit_boundary[ipar] = false;
                        if (m_logger != NULL) {
                            *m_logger << "... parameter \"";
                            *m_logger << pars->par(ipar)->name();
                            *m_logger << "\" does not drive optimization step anymore.";
                            *m_logger << std::endl;
                        }
                    }
                }
            }

        } // endfor: looped over all parameters

        // Signal if a parameter is driving the optimization step
        if (ipar_bnd != -1) {
            m_hit_boundary[ipar_bnd] = true;
            if (m_logger != NULL) {
                *m_logger << "... parameter \"";
                *m_logger << pars->par(ipar_bnd)->name();
                *m_logger << "\" drives optimization step (step=";
                *m_logger << step << ")" << std::endl;
            }
        }

    } // endif: automatic step size adjustment requested

    // Return step size
    return step;
}


/***********************************************************************//**
 * @brief Compute parameter uncertainties
 *
 * @param[in] fct Optimizer function.
 * @param[in] pars Function parameters.
 *
 * Compute parameter uncertainties from the diagonal elements of the
 * covariance matrix.
 ***************************************************************************/
void GOptimizerLM::errors(GOptimizerFunction* fct, GOptimizerPars* pars)
{
    // Get number of parameters
    int npars = pars->npars();

    // Perform final parameter evaluation
    fct->eval(*pars);

    // Fetch sparse matrix pointer. We have to do this after the eval()
    // method since eval() will allocate new memory for the covariance
    // matrix!
    GSparseMatrix* covar = fct->covar();

    // Save best fitting value
    m_value = *(fct->value());

    // Save covariance matrix
    GSparseMatrix save_covar = GSparseMatrix(*covar);

    // Signal no diagonal element loading
    bool diag_loaded = false;

    // Loop over error computation (maximum 2 turns)
    for (int i = 0; i < 2; ++i) {

        // Solve: covar * X = unit
        try {
            covar->cholesky_decompose(1);
            GVector unit(npars);
            for (int ipar = 0; ipar < npars; ++ipar) {
                unit(ipar) = 1.0;
                GVector x  = covar->cholesky_solver(unit,1);
                if (x(ipar) >= 0.0)
                    pars->par(ipar)->error(sqrt(x(ipar)));
                else {
                    pars->par(ipar)->error(0.0);
                    m_status = G_LM_BAD_ERRORS;
                }
                unit(ipar) = 0.0;
            }
        }
        catch (GException::matrix_zero &e) {
            m_status = G_LM_SINGULAR;
            if (m_logger != NULL) {
                *m_logger << "GOptimizerLM::terminate: "
                          << "All curvature matrix elements are zero."
                          << std::endl;
            }
            break;
        }
        catch (GException::matrix_not_pos_definite &e) {

            // Load diagonal if this has not yet been tried
            if (!diag_loaded) {

                // Flag errors as inaccurate
                m_status = G_LM_BAD_ERRORS;
                if (m_logger != NULL) {
                    *m_logger << "Non-Positive definite curvature matrix encountered."
                              << std::endl;
                    *m_logger << "Load diagonal elements with 1e-10."
                              << " Fit errors may be inaccurate."
                              << std::endl;
                }

                // Try now with diagonal loaded matrix
                *covar = save_covar;
                for (int ipar = 0; ipar < npars; ++ipar)
                    (*covar)(ipar,ipar) += 1.0e-10;

                // Signal loading
                diag_loaded = true;

                // Try again
                continue;

            } // endif: diagonal has not yet been loaded

            // ... otherwise signal an error
            else {
                m_status = G_LM_NOT_POSTIVE_DEFINITE;
                if (m_logger != NULL) {
                    *m_logger << "Non-Positive definite curvature matrix encountered,"
                              << " even after diagonal loading." << std::endl;
                }
                break;
            }
        }
        catch (std::exception &e) {
            throw;
        }

        // If no error occured then break now
        break;

    } // endfor: looped over error computation

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put optimizer in output stream
 *
 * @param[in] os Output stream into which the optimizer will be dumped.
 * @param[in] opt Object to be dumped.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GOptimizerLM& opt)
{
    // Put optimizer in stream
    os << "=== GOptimizerLM ===" << std::endl;
    os << " Optimized function value ..: " << opt.m_value << std::endl;
    os << " Absolute precision ........: " << opt.m_eps << std::endl;
    os << " Optimization status .......: " << opt.m_status << std::endl;
    os << " Number of iterations ......: " << opt.m_iter << std::endl;
    os << " Lambda ....................: " << opt.m_lambda;

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Write optimizer into logger
 *
 * @param[in] log Logger
 * @param[in] opt Optimizer to be writted into logger.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GOptimizerLM& opt)
{
    // Write optimizer into logger
    log << "=== GOptimizerLM ===" << std::endl;
    log << " Optimized function value ..: " << opt.m_value << std::endl;
    log << " Absolute precision ........: " << opt.m_eps << std::endl;
    log << " Optimization status .......: " << opt.m_status << std::endl;
    log << " Number of iterations ......: " << opt.m_iter << std::endl;
    log << " Lambda ....................: " << opt.m_lambda;

    // Return logger
    return log;
}
