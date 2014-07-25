/***************************************************************************
 *             GOptimizerLM.cpp - Levenberg Marquardt optimizer            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2014 by Juergen Knoedlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GOptimizerLM.cpp
 * @brief Levenberg-Marquardt optimizer class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GOptimizerLM.hpp"
#include "GTools.hpp"
#include "GException.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_OPT              //!< Define to debug optimize() method
//#define G_DEBUG_ITER             //!< Define to debug iteration() method
//#define G_DEBUG_SHOW_GRAD_COVAR  //!< Define to show grad and curvature


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
GOptimizerLM::~GOptimizerLM(void)
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


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear object
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GOptimizerLM::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GOptimizer::free_members();

    // Initialise members
    this->GOptimizer::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone object
***************************************************************************/
GOptimizerLM* GOptimizerLM::clone(void) const
{
    return new GOptimizerLM(*this);
}


/***********************************************************************//**
 * @brief Optimize function parameters
 *
 * @param[in] fct Optimization function.
 * @param[in] pars Function parameters.
 ***************************************************************************/
void GOptimizerLM::optimize(GOptimizerFunction& fct, GOptimizerPars& pars)
{
    // Get number of parameters. Continue only if there are free parameters
    m_npars = pars.size();
    m_nfree = pars.nfree();
    if (m_nfree > 0) {

        // Initialise optimization parameters
        m_lambda  = m_lambda_start;
        m_status  = G_LM_CONVERGED;

        // Initialise bookkeeping arrays
        m_hit_boundary.clear();
        m_hit_minimum.clear();
        m_hit_maximum.clear();
        m_par_freeze.clear();
        m_par_remove.clear();
        m_hit_boundary.reserve(m_npars);
        m_hit_minimum.reserve(m_npars);
        m_hit_maximum.reserve(m_npars);
        m_par_freeze.reserve(m_npars);
        m_par_remove.reserve(m_npars);
        for (int i = 0; i < m_npars; ++i) {
            m_hit_boundary.push_back(false);
            m_hit_minimum.push_back(0);
            m_hit_maximum.push_back(0);
            m_par_freeze.push_back(false);
            m_par_remove.push_back(false);
        }

        // Initial function evaluation
        fct.eval(pars);

        // If a free parameter has a zero diagonal element in the curvature
        // matrix then remove this parameter definitely from the fit as it
        // otherwise will block the fit. The problem appears in the unbinned
        // fitting where parameter gradients may be zero (due to the truncation
        // of the PSF), but the Npred gradient is not zero. In principle we
        // could use the Npred gradient for fitting (I guess), but I still
        // have to figure out how ... (the diagonal loading was not so
        // successful as it faked early convergence)
        for (int ipar = 0; ipar < m_npars; ++ipar) {
            if (pars[ipar]->is_free()) {
                if ((*fct.curvature())(ipar,ipar) == 0.0) {
                    if (m_logger != NULL) {
                        *m_logger << "  Parameter \"" << pars[ipar]->name();
                        *m_logger << "\" has zero curvature.";
                        *m_logger << " Fix parameter." << std::endl;
                    }
                    m_par_remove[ipar] = true;
                    pars[ipar]->fix();
                }
            }
        }

        // Save function value
        m_value = fct.value();

        // Save initial statistics and lambda values
        double value_old  = m_value;
        double lambda_old = m_lambda;
        int    lambda_inc = 0;

        // Optionally write initial iteration into logger
        if (m_logger != NULL) {
            (*m_logger)(">Iteration %3d: -logL=%.3f, Lambda=%.1e",
                        0, m_value, m_lambda);   
        }
        #if defined(G_DEBUG_OPT)
        std::cout << "Initial iteration: func=" << m_value << ", Lambda="
                  << m_lambda << std::endl;
        #endif
        #if defined(G_DEBUG_SHOW_GRAD_COVAR)
        if (m_logger != NULL) {
            *m_logger << *fct.gradient() << std::endl;
            *m_logger << *fct.curvature() << std::endl;
        }
        #endif

        // Initialise iteration flags
        bool check_for_freeze = true;

        // Iterative fitting
        for (m_iter = 1; m_iter <= m_max_iter; ++m_iter) {

            // Perform one iteration
            iteration(fct, pars);

            // Compute function improvement (>0 means decrease)
            //double delta = value_old - m_value;

            // Determine maximum (scaled) gradient
            double grad_max  = 0.0;
            int    grad_imax = -1;
            for (int ipar = 0; ipar < m_npars; ++ipar) {
                if (pars[ipar]->is_free()) {
                    double grad = pars[ipar]->factor_gradient();
                    if (std::abs(grad) > std::abs(grad_max)) {
                        grad_max  = grad;
                        grad_imax = ipar;
                    }
                    if (grad == 0.0) {
                        if (m_logger != NULL) {
                            *m_logger << "Parameter " << ipar;
                            *m_logger << " (" << pars[ipar]->name() << ")";
                            *m_logger << " has a zero gradient." << std::endl;
                        }
                    }
                }
            }

            // Optionally write iteration results into logger
            if (m_logger != NULL) {
                std::string stalled = "";
                std::string status  = "";
                if (m_lambda > lambda_old) {
                    status  = " ";
                    stalled = " (stalled)";
                }
                else {
                    status = ">";
                }
                std::string parname = "";
                if (grad_imax != -1) {
                    parname = " [" + pars[grad_imax]->name() + ":" +
                              gammalib::str(grad_imax) + "]";
                }
                (*m_logger)("%sIteration %3d: -logL=%.3f, Lambda=%.1e,"
                            " delta=%.3f, max(|grad|)=%f%s%s",
                            status.c_str(), m_iter, m_value, lambda_old,
                            m_delta, grad_max,
                            parname.c_str(), stalled.c_str());   
            }
            #if defined(G_DEBUG_OPT)
            std::cout << "Iteration " << m_iter << ": func=" 
                      << m_value << ", Lambda=" << lambda_old
                      << ", delta=" << m_delta << std::endl;
            #endif
            #if defined(G_DEBUG_SHOW_GRAD_COVAR)
            if (m_logger != NULL) {
                *m_logger << *fct.gradient() << std::endl;
                *m_logger << *fct.curvature() << std::endl;
            }
            #endif

            // Reset lambda increment if we had success
            if (m_lambda < lambda_old) {
                lambda_inc = 0;
            }

            // If function increased while lambda did not increase then stop
            // iterations
            //if ((m_lambda <= lambda_old) && (delta < 0.0))
            //    break;

            // Stop if convergence was reached. Before stopping, check
            // if some parameters were frozen, and if this was the case,
            // free them now and continue. We do this only once, i.e.
            // the next time a parameter is frozen and convergence is
            // reached we really stop.
            if ((m_lambda <= lambda_old) && (m_delta < m_eps)) {

                // Check for frozen parameters, and if some exist, free
                // them and start over
                if (check_for_freeze) {

                    // Signal that we won't check frozen parameters
                    // anymore
                    check_for_freeze = false;

                    // Free frozen parameters and determine how many
                    // have been frozen
                    int nfrozen = 0;
                    for (int ipar = 0; ipar < m_npars; ++ipar) {
                        if (m_par_freeze[ipar]) {
                            nfrozen++;
                            pars[ipar]->free();
                            if (m_logger != NULL) {
                                *m_logger << "  Free parameter \""
                                          << pars[ipar]->name()
                                          << "\" after convergence was"
                                             " reached with frozen"
                                             " parameter." << std::endl;
                            }
                        }
                    }

                    // If there were frozen parameters then start over
                    // again (initialise optimizer)
                    if (nfrozen > 0) {
                        m_lambda   = m_lambda_start;
                        value_old  = m_value;
                        lambda_old = m_lambda;
                        lambda_inc = 0;
                        continue;
                    }
                }

                // ... otherwise, convergence was reached and we can stop
                // now
                break;
                
            } // endif: convergence check

            // Monitor the number of subsequent increases of lambda and
            // stop if the number of increases exceeds threshold
            lambda_inc = (m_lambda > lambda_old) ? lambda_inc + 1 : 0;
            if (lambda_inc >= m_max_stall) {
                m_status = G_LM_STALLED;
                break;
            }

            // Bookkeeping of actual result (we always store the last
            // lambda to detect turn arounds in the lambda tendency; however
            // we always keep the best function value)
            lambda_old = m_lambda;
            if (m_delta > 0.0) {
                value_old = m_value;
            }

        } // endfor: iterations

        // Compute parameter uncertainties
        errors(fct, pars);
        
        // Free now all temporarily frozen parameters so that the resulting
        // model has the same attributes as the initial model
        for (int ipar = 0; ipar < m_npars; ++ipar) {
            if (m_par_freeze[ipar] || m_par_remove[ipar]) {
                pars[ipar]->free();
                if (m_logger != NULL) {
                    *m_logger << "  Free parameter \""
                              << pars[ipar]->name()
                              << "\" after convergence was"
                                 " reached with frozen"
                                 " parameter." << std::endl;
                }
            }
        }

    } // endif: there were free parameters to fit

    // ... otherwise just execute final step
    else {
        errors(fct, pars);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print optimizer information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing optimizer information.
 ***************************************************************************/
std::string GOptimizerLM::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GOptimizerLM ===");

        // Append information
        result.append("\n"+gammalib::parformat("Optimized function value"));
        result.append(gammalib::str(m_value, 3));
        result.append("\n"+gammalib::parformat("Absolute precision"));
        result.append(gammalib::str(m_eps));

        // Append status
        result.append("\n"+gammalib::parformat("Optimization status"));
        switch (m_status) {
        case G_LM_CONVERGED:
            result.append("converged");
            break;
        case G_LM_STALLED:
            result.append("stalled");
            break;
        case G_LM_SINGULAR:
            result.append("singular curvature matrix encountered");
            break;
        case G_LM_NOT_POSTIVE_DEFINITE:
            result.append("curvature matrix not positive definite");
            break;
        case G_LM_BAD_ERRORS:
            result.append("errors are inaccurate");
            break;
        default:
            result.append("unknown");
            break;
        }

        // Append further information
        result.append("\n"+gammalib::parformat("Number of parameters"));
        result.append(gammalib::str(m_npars));
        result.append("\n"+gammalib::parformat("Number of free parameters"));
        result.append(gammalib::str(m_nfree));
        result.append("\n"+gammalib::parformat("Number of iterations"));
        result.append(gammalib::str(m_iter));
        result.append("\n"+gammalib::parformat("Lambda"));
        result.append(gammalib::str(m_lambda));

    } // endif: chatter was not silent

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
void GOptimizerLM::init_members(void)
{
    // Initialise optimizer parameters
    m_npars        = 0;
    m_nfree        = 0;
    m_lambda_start = 1.0e-3;
    m_lambda_inc   = 10.0;
    m_lambda_dec   = 0.1;
    m_eps          = 1.0e-6;
    m_max_iter     = 1000;
    m_max_stall    = 10;
    m_max_hit      = 3; //!< Maximum successive boundary hits before freeze
    m_step_adjust  = true;

    // Initialise bookkeeping arrays
    m_hit_boundary.clear();
    m_hit_minimum.clear();
    m_hit_maximum.clear();
    m_par_freeze.clear();
    m_par_remove.clear();

    // Initialise optimizer values
    m_lambda = m_lambda_start;
    m_value  = 0.0;
    m_delta  = 0.0;
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
    m_hit_boundary = opt.m_hit_boundary;
    m_hit_minimum  = opt.m_hit_minimum;
    m_hit_maximum  = opt.m_hit_maximum;
    m_par_freeze   = opt.m_par_freeze;
    m_par_remove   = opt.m_par_remove;
    m_lambda       = opt.m_lambda;
    m_value        = opt.m_value;
    m_delta        = opt.m_delta;
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
void GOptimizerLM::iteration(GOptimizerFunction& fct, GOptimizerPars& pars)
{
    // Debug option: dump header
    #if defined(G_DEBUG_ITER)
    std::cout << "GOptimizerLM::iteration: enter" << std::endl;
    #endif

    // Single loop for common exit point
    do {

        // Initialise iteration parameters
        GVector*       grad      = fct.gradient();
        GMatrixSparse* curvature = fct.curvature();

        // Save function value, gradient and curvature matrix
        double         save_value     = m_value;
        GVector        save_grad      = GVector(*grad);
        GMatrixSparse  save_curvature = GMatrixSparse(*curvature);

        // Save parameter values in vector
        GVector save_pars(m_npars);
        for (int ipar = 0; ipar < m_npars; ++ipar) {
            save_pars[ipar] = pars[ipar]->factor_value();
        }

        // Setup matrix and vector for curvature computation
        for (int ipar = 0; ipar < m_npars; ++ipar) {
            (*curvature)(ipar,ipar) *= (1.0 + m_lambda);
            (*grad)[ipar]            = -(*grad)[ipar];
        }
        
        // Debug option: dump gradient and curvature matrix
        #if defined(G_DEBUG_ITER)
        std::cout << "Gradient : ";
        for (int ipar = 0; ipar < m_npars; ++ipar) {
            std::cout << (*grad)[ipar] << " ";
        }
        std::cout << std::endl;
        std::cout << "Curvature: ";
        for (int ipar = 0; ipar < m_npars; ++ipar) {
            for (int jpar = 0; jpar < m_npars; ++jpar) {
                std::cout << (*curvature)(ipar,jpar) << " ";
            }
        }
        std::cout << std::endl;
        #endif

        // Solve: curvature * X = grad. Handle matrix problems
        try {
            *grad = curvature->solve(*grad);
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
        double step = step_size(*grad, pars);

        // Debug option: dump step size
        #if defined(G_DEBUG_ITER)
        std::cout << "Step size = " << step << std::endl;
        #endif

        // Derive new parameter vector
        for (int ipar = 0; ipar < m_npars; ++ipar) {

            // Consider only free parameters
            if (pars[ipar]->is_free()) {

                // Get actual parameter value and limits
                double p     = pars[ipar]->factor_value();
                double p_min = pars[ipar]->factor_min();
                double p_max = pars[ipar]->factor_max();

                // Compute new parameter value
                p += (*grad)[ipar] * step;

                // Debug option: dump new parameter value
                #if defined(G_DEBUG_ITER)
                std::cout << "Trial factor value of ";
                std::cout << pars[ipar]->name() << " = ";
                std::cout << p << std::endl;
                #endif

                // Constrain parameter to within the valid range
                if (pars[ipar]->has_min() && p < p_min) {
                    if (m_hit_minimum[ipar] >= m_max_hit) {
                        if (m_logger != NULL) {
                            *m_logger << "  Parameter \"" << pars[ipar]->name();
                            *m_logger << "\" hits minimum " << p_min << " more than ";
                            *m_logger << m_max_hit << " times.";
                            *m_logger << " Fix parameter at minimum for now." << std::endl;
                        }
                        m_par_freeze[ipar] = true;
                        pars[ipar]->fix();
                    }
                    else {
                        if (m_logger != NULL) {
                            *m_logger << "  Parameter \"" << pars[ipar]->name();
                            *m_logger << "\" hits minimum: ";
                            *m_logger << p << " < " << p_min;
                            *m_logger << " (" << m_hit_minimum[ipar]+1 << ")" << std::endl;
                        }
                    }
                    m_hit_minimum[ipar]++;
                    p = p_min;
                }
                else if (pars[ipar]->has_max() && p > p_max) {
                    if (m_hit_maximum[ipar] >= m_max_hit) {
                        if (m_logger != NULL) {
                            *m_logger << "  Parameter \"" << pars[ipar]->name();
                            *m_logger << "\" hits maximum " << p_max << " more than ";
                            *m_logger << m_max_hit << " times.";
                            *m_logger << " Fix parameter at maximum for now." << std::endl;
                        }
                        m_par_freeze[ipar] = true;
                        pars[ipar]->fix();
                    }
                    else {
                        if (m_logger != NULL) {
                            *m_logger << "  Parameter \"" << pars[ipar]->name();
                            *m_logger << "\" hits maximum: ";
                            *m_logger << p << " > " << p_max;
                            *m_logger << " (" << m_hit_maximum[ipar]+1 << ")" << std::endl;
                        }
                    }
                    m_hit_maximum[ipar]++;
                    p = p_max;
                }
                else {
                    m_hit_minimum[ipar] = 0;
                    m_hit_maximum[ipar] = 0;
                }

                // Set new parameter value
                pars[ipar]->factor_value(p);

                // Debug option: dump new parameter value
                #if defined(G_DEBUG_ITER)
                std::cout << "New value of ";
                std::cout << pars[ipar]->name() << " = ";
                std::cout << pars[ipar]->value() << std::endl;
                #endif

            } // endif: Parameter was free

        } // endfor: computed new parameter vector

        // Evaluate function at new parameters
        fct.eval(pars);

        // If a free parameter has a zero diagonal element in the curvature
        // matrix then remove this parameter definitely from the fit as it
        // otherwise will block the fit. The problem appears in the unbinned
        // fitting where parameter gradients may be zero (due to the truncation
        // of the PSF), but the Npred gradient is not zero. In principle we
        // could use the Npred gradient for fitting (I guess), but I still
        // have to figure out how ... (the diagonal loading was not so
        // successful as it faked early convergence)
        for (int ipar = 0; ipar < m_npars; ++ipar) {
            if (pars[ipar]->is_free()) {
                if ((*fct.curvature())(ipar,ipar) == 0.0) {
                    if (m_logger != NULL) {
                        *m_logger << "  Parameter \"" << pars[ipar]->name();
                        *m_logger << "\" has zero curvature.";
                        *m_logger << " Fix parameter." << std::endl;
                    }
                    m_par_remove[ipar] = true;
                    pars[ipar]->fix();
                }
            }
        }

        // Fetch new pointers since eval will allocate new memory
        grad      = fct.gradient();
        curvature = fct.curvature();

        // Retrieve new function value
        m_value = fct.value();

        // Debug option: dump new function value
        #if defined(G_DEBUG_ITER)
        std::cout << "Function value " << m_value;
        std::cout << " (was " << save_value << ")" << std::endl;
        #endif

        // Determine how many parameters have changed
        int par_change = 0;
        for (int ipar = 0; ipar < m_npars; ++ipar) {
            if (pars[ipar]->factor_value() != save_pars[ipar]) {
                par_change++;
            }
        }

        // If the function has decreased then accept the new solution
        // and decrease lambda ...
        m_delta = save_value - m_value;
        if (m_delta > 0.0) {
            m_lambda *= m_lambda_dec;
        }

        // ... if function is identical then accept new solution. If the
        // parameters have changed then increase lambda, otherwise decrease
        // lambda
        else if (m_delta == 0.0) {
            if (par_change) {
                m_lambda *= m_lambda_inc;
            }
            else {
                m_lambda *= m_lambda_dec;
            }
        }

        // ... if function worsened slightly then accept new solution and
        // increase lambda
        else if (m_delta > -1.0e-6) {
            m_lambda *= m_lambda_inc;
        }

        // ... otherwise, if the statistics did not improve then use old
        // parameters and increase lamdba. Restore also the best statistics
        // value that was reached so far, the gradient vector and the curve
        // matrix.
        else {
            m_lambda  *= m_lambda_inc;
            m_value    = save_value;
            *grad      = save_grad;
            *curvature = save_curvature;
            for (int ipar = 0; ipar < m_npars; ++ipar) {
                pars[ipar]->factor_value(save_pars[ipar]);
            }
        }

    } while (0); // endwhile: main loop

    // Debug option: dump trailer
    #if defined(G_DEBUG_ITER)
    std::cout << "GOptimizerLM::iteration: exit" << std::endl;
    #endif

    // Return
    return;

}


/***********************************************************************//**
 * @brief Compute parameter uncertainties
 *
 * @param[in] fct Optimizer function.
 * @param[in] pars Function parameters.
 *
 * Compute parameter uncertainties from the diagonal elements of the
 * curvature matrix.
 ***************************************************************************/
void GOptimizerLM::errors(GOptimizerFunction& fct, GOptimizerPars& pars)
{
    // Get number of parameters
    int npars = pars.size();

    // Perform final parameter evaluation
    fct.eval(pars);

    // Fetch sparse matrix pointer. We have to do this after the eval()
    // method since eval() will allocate new memory for the curvature
    // matrix!
    GMatrixSparse* curvature = fct.curvature();

    // Save best fitting value
    m_value = fct.value();

    // Save curvature matrix
    GMatrixSparse save_curvature = GMatrixSparse(*curvature);

    // Signal no diagonal element loading
    bool diag_loaded = false;

    // Loop over error computation (maximum 2 turns)
    for (int i = 0; i < 2; ++i) {

        // Solve: curvature * X = unit
        try {
            GMatrixSparse decomposition = curvature->cholesky_decompose(true);
            GVector unit(npars);
            for (int ipar = 0; ipar < npars; ++ipar) {
                unit[ipar] = 1.0;
                GVector x  = decomposition.cholesky_solver(unit, true);
                if (x[ipar] >= 0.0) {
                    pars[ipar]->factor_error(sqrt(x[ipar]));
                }
                else {
                    pars[ipar]->factor_error(0.0);
                    m_status = G_LM_BAD_ERRORS;
                }
                unit[ipar] = 0.0;
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
                *curvature = save_curvature;
                for (int ipar = 0; ipar < npars; ++ipar) {
                    (*curvature)(ipar,ipar) += 1.0e-10;
                }

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
double GOptimizerLM::step_size(const GVector& grad, const GOptimizerPars& pars)
{
    // Initialise step size
    double step = 1.0;

    // Check if we should reduce the step size
    if (m_step_adjust) {

        // Initialise the parameter index that constrains most the fit
        int ipar_bnd = -1;

        // Loop over all parameters
        for (int ipar = 0; ipar < pars.size(); ++ipar) {

            // Get parameter attributes
            double p     = pars[ipar]->factor_value();
            double p_min = pars[ipar]->factor_min();
            double p_max = pars[ipar]->factor_max();
            double delta = grad[ipar];

            // Check if a parameter minimum requires a reduced step size
            if (pars[ipar]->has_min()) {
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
                            *m_logger << "  Parameter \"";
                            *m_logger << pars[ipar]->name();
                            *m_logger << "\" does not drive optimization step anymore.";
                            *m_logger << std::endl;
                        }
                    }
                }
            }

            // Check if a parameter maximum requires a reduced step size
            if (pars[ipar]->has_max()) {
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
                            *m_logger << "  Parameter \"";
                            *m_logger << pars[ipar]->name();
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
                *m_logger << "  Parameter \"";
                *m_logger << pars[ipar_bnd]->name();
                *m_logger << "\" drives optimization step (step=";
                *m_logger << step << ")" << std::endl;
            }
        }

    } // endif: automatic step size adjustment requested

    // Return step size
    return step;
}
