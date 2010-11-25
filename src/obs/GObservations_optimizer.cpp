/***************************************************************************
 *  GObservations_optimizer.cpp  -  Optimizer class of observations class  *
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
 * @file GObservations_optimizer.cpp
 * @brief GObservations::optimizer class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GObservations.hpp"
#include "GEventList.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
#define G_EVAL_TIMING 0 //!< Perform optimizer timing (0=no, 1=yes)
#define G_EVAL_DEBUG  0 //!< Perform optimizer debugging (0=no, 1=yes)

/* __ Prototypes _________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
***************************************************************************/
GObservations::optimizer::optimizer(void) : GOptimizerFunction()
{
    // Initialise iterator
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor based on specific object
 *
 * @param[in] obs Observations container to use for optimizer.
 ***************************************************************************/
GObservations::optimizer::optimizer(GObservations *obs) : GOptimizerFunction()
{
    // Initialise iterator
    init_members();

    // Set object
    m_this = obs;

    // Return
    return;
}



/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] obs Instance which should be used for construction
 ***************************************************************************/
GObservations::optimizer::optimizer(const optimizer& fct) : GOptimizerFunction(fct)
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
GObservations::optimizer::~optimizer()
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
 * @param[in] obs Instance to be assigned
 ***************************************************************************/
GObservations::optimizer& GObservations::optimizer::operator= (const optimizer& fct)
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
 =                               Public methods                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Evaluate log-likelihood function
 *
 * @param[in] pars Optimizer parameters.
 *
 * This method evaluates the log-likelihood function for parameter
 * optimisation. It handles both binned and unbinned data. 
 * For binned data the function to optimize is given by
 * \f$L=-\sum_i n_i \log e_i - e_i\f$
 * where the sum is taken over all data space bins, \f$n_i\f$ is the
 * observed number of counts and \f$e_i\f$ is the model.
 * For unbinned data the function to optimize is given by
 * \f$L=-\sum_i \log e_i + {\rm Npred}\f$
 * where the sum is taken over all events and \f${\rm Npred}\f$ is the
 * total number of events that is predicted by the model.
 * Note that binned and unbinned observations may be combined.
 ***************************************************************************/
/*
void GObservations::optimizer::eval(const GOptimizerPars& pars) 
{
    // Timing measurement
    #if G_EVAL_TIMING
    clock_t t_start = clock();
    #endif

    // Allocate pointers for temporary memory
    int*    inx    = NULL;
    double* values = NULL;

    // Single loop for common exit point
    do {
        // Get number of parameters
        int npars = pars.npars();

        // Fall through if we have no free parameters
        if (npars < 1)
            continue;

        // Free old memory
        if (m_gradient != NULL) delete m_gradient;
        if (m_covar    != NULL) delete m_covar;

        // Initialise value, gradient vector and curvature matrix
        m_value    = 0.0;
        m_npred    = 0.0;
        m_gradient = new GVector(npars);
        m_covar    = new GSparseMatrix(npars,npars);

        // Allocate some working arrays
        GVector grad(npars);
        m_covar->stack_init(npars,10000);
        inx    = new int[npars];
        values = new double[npars];

        // Collect predicted number of events for unbinned observations
        for (int i = 0; i < m_this->m_num; ++i) {
            if (m_this->m_obs[i]->events()->islist()) {
                m_npred += m_this->m_obs[i]->npred((GModels&)pars, &grad);
                for (int k = 0; k < npars; ++k)
                    (*m_gradient)(k) += grad(k);
                #if G_EVAL_DEBUG
                std::cout << "Npred=" << m_npred << " Grad="
                          << *m_gradient << std::endl;
                #endif
            }
        }
        m_value += m_npred;

        // Iterate over all data bins
        GObservations::iterator end = m_this->end();
        for (GObservations::iterator bin = m_this->begin(); bin != end; ++bin) {

            // Get number of counts in bin
            double data = bin->counts();

            // Get model and derivative
            double model = bin->model((GModels&)pars, &grad);

            // Skip bin if model is too small (avoids -Inf or NaN gradients)
            if (model <= m_minmod)
                continue;

            // Create index array of non-zero derivatives
            int ndev = 0;
            for (int i = 0; i < npars; ++i) {
                if (grad(i) != 0.0 && !std::isinf(grad(i))) {
                    inx[ndev] = i;
                    ndev++;
                }
            }

            // Case A: binned analysis
            if (bin->isbin()) {

                // Update Poissonian statistics (excluding factorial
                // term for faster computation)
                m_value -= data * log(model) - model;

                // Skip bin now if there are no non-zero derivatives
                if (ndev < 1)
                    continue;

                // Update gradient vector and curvature matrix. To avoid
                // unneccessary computations we distinguish the case where
                // data>0 and data=0. The second case requires much less
                // computation since it does not contribute to the covariance
                // matrix ...
                if (data > 0.0) {

                    // Pre computation
                    double fb = data / model;
                    double fc = (1.0 - fb);
                    double fa = fb / model;

                    // Loop over columns
                    for (int jdev = 0; jdev < ndev; ++jdev) {

                        // Initialise computation
                        register int jpar    = inx[jdev];
                        double       g       = grad(jpar);
                        double       fa_i    = fa * g;

                        // Update gradient
                        (*m_gradient)(jpar) += fc * g;

                        // Loop over rows
                        register int* ipar = inx;
                        for (register int idev = 0; idev < ndev; ++idev, ++ipar)
                            values[idev] = fa_i * grad(*ipar);

                        // Add column to matrix
                        m_covar->add_col(values, inx, ndev, jpar);
                    } // endfor: looped over columns
                } // endif: data was > 0

                // ... handle now data=0
                else {
                    register int* ipar = inx;
                    for (register int idev = 0; idev < ndev; ++idev, ++ipar)
                        (*m_gradient)(*ipar) += grad(*ipar);
                }

            } // endif: analysis was binned

            // Case B: unbinned analysis
            else {

                // Update Poissonian statistics (excluding factorial term
                // for faster computation)
                m_value -= log(model);

                // Skip bin now if there are no non-zero derivatives
                if (ndev < 1)
                    continue;

                // Update gradient vector and curvature matrix.
                double fb = 1.0 / model;
                double fa = fb / model;
                for (int jdev = 0; jdev < ndev; ++jdev) {

                    // Initialise computation
                    register int jpar    = inx[jdev];
                    double       g       = grad(jpar);
                    double       fa_i    = fa * g;

                    // Update gradient.
                    (*m_gradient)(jpar) -= fb * g;

                    // Loop over rows
                    register int* ipar = inx;
                    for (register int idev = 0; idev < ndev; ++idev, ++ipar)
                        values[idev] = fa_i * grad(*ipar);

                    // Add column to matrix
                    m_covar->add_col(values, inx, ndev, jpar);

                } // endfor: looped over columns

            } // endelse: analysis was unbinned

        } // endfor: iterated over all data bins

        // Release stack
        m_covar->stack_destroy();

    } while(0); // endwhile: main loop

    // Free temporary memory
    if (values != NULL) delete [] values;
    if (inx    != NULL) delete [] inx;

    // Optionally dump gradient and covariance matrix
    #if G_EVAL_DEBUG
    std::cout << *m_gradient << std::endl;
    std::cout << *m_covar << std::endl;
    #endif

    // Timing measurement
    #if G_EVAL_TIMING
    double t_elapse = (double)(clock() - t_start) / (double)CLOCKS_PER_SEC;
    std::cout << "GObservations::optimizer::eval: CPU usage = "
              << t_elapse << " sec" << std::endl;
    #endif

    // Return
    return;
}
*/

/***********************************************************************//**
 * @brief Evaluate log-likelihood function
 *
 * @param[in] pars Optimizer parameters.
 *
 * This method evaluates the log-likelihood function for parameter
 * optimisation. It handles both binned and unbinned data. 
 * For binned data the function to optimize is given by
 * \f$L=-\sum_i n_i \log e_i - e_i\f$
 * where the sum is taken over all data space bins, \f$n_i\f$ is the
 * observed number of counts and \f$e_i\f$ is the model.
 * For unbinned data the function to optimize is given by
 * \f$L=-\sum_i \log e_i + {\rm Npred}\f$
 * where the sum is taken over all events and \f${\rm Npred}\f$ is the
 * total number of events that is predicted by the model.
 * Note that binned and unbinned observations may be combined.
 ***************************************************************************/
void GObservations::optimizer::eval(const GOptimizerPars& pars) 
{
    // Timing measurement
    #if G_EVAL_TIMING
    clock_t t_start = clock();
    #endif

    // Single loop for common exit point
    do {
        // Get number of parameters
        int npars = pars.npars();

        // Fall through if we have no free parameters
        if (npars < 1)
            continue;

        // Free old memory
        if (m_gradient != NULL) delete m_gradient;
        if (m_covar    != NULL) delete m_covar;
        if (m_wrk_grad != NULL) delete m_wrk_grad;

        // Initialise value, gradient vector and curvature matrix
        m_value    = 0.0;
        m_npred    = 0.0;
        m_gradient = new GVector(npars);
        m_covar    = new GSparseMatrix(npars,npars);
        m_wrk_grad = new GVector(npars);
        m_covar->stack_init(npars,10000);

        // Loop over all observations
        for (int i = 0; i < m_this->m_num; ++i) {

            // Unbinned analysis
            if (m_this->m_obs[i]->events()->islist()) {

                // Update the Npred value
                m_npred     += m_this->m_obs[i]->npred((GModels&)pars, m_wrk_grad);
                *m_gradient += *m_wrk_grad;
                #if G_EVAL_DEBUG
                std::cout << "Npred=" << m_npred << " Grad="
                          << *m_gradient << std::endl;
                #endif

                // Update the log-likelihood
                poisson_unbinned(*(m_this->m_obs[i]), pars);
                
            } // endif: unbinned analysis

            // ... or binned analysis
            else {

                // Update the log-likelihood
                poisson_binned(*(m_this->m_obs[i]), pars);
            
            } // endelse: binned analysis

        }

        // Add the Npred value to the log-likelihood
        m_value += m_npred;

        // Release stack
        m_covar->stack_destroy();

    } while(0); // endwhile: main loop

    // Optionally dump gradient and covariance matrix
    #if G_EVAL_DEBUG
    std::cout << *m_gradient << std::endl;
    std::cout << *m_covar << std::endl;
    #endif

    // Timing measurement
    #if G_EVAL_TIMING
    double t_elapse = (double)(clock() - t_start) / (double)CLOCKS_PER_SEC;
    std::cout << "GObservations::optimizer::eval: CPU usage = "
              << t_elapse << " sec" << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Evaluate log-likelihood function for Poisson statistics and
 * unbinned analysis
 *
 * @param[in] obs Observation.
 * @param[in] pars Optimizer parameters.
 *
 * Evaluate log-likelihood function, gradient and curvature matrix for one
 * observation using unbinned analysis and Poisson statistics.
 *
 * @todo We still use the pointing information from the event. This has to
 * be removed later ...
 ***************************************************************************/
void GObservations::optimizer::poisson_unbinned(const GObservation& obs,
                                                const GOptimizerPars& pars) 
{
    // Timing measurement
    #if G_EVAL_TIMING
    clock_t t_start = clock();
    #endif

    // Get number of parameters
    int npars = pars.npars();

    // Allocate some working arrays
    int*    inx    = new int[npars];
    double* values = new double[npars];

    // Iterate over all events
    for (int i = 0; i < obs.events()->size(); ++i) {

        // Get pointer to event
        GEvent* event = obs.events()->pointer(i);

        // Get model and derivative
        double model = obs.model((GModels&)pars, *(event->dir()),
                                 *(event->energy()), *(event->time()),
                                 *(event->pnt()), m_wrk_grad);

        // Skip bin if model is too small (avoids -Inf or NaN gradients)
        if (model <= m_minmod)
            continue;

        // Create index array of non-zero derivatives
        int ndev = 0;
        for (int i = 0; i < npars; ++i) {
            if ((*m_wrk_grad)(i) != 0.0 && !std::isinf((*m_wrk_grad)(i))) {
                inx[ndev] = i;
                ndev++;
            }
        }

        // Update Poissonian statistics (excluding factorial term for faster
        // computation)
        m_value -= log(model);

        // Skip bin now if there are no non-zero derivatives
        if (ndev < 1)
            continue;

        // Update gradient vector and curvature matrix.
        double fb = 1.0 / model;
        double fa = fb / model;
        for (int jdev = 0; jdev < ndev; ++jdev) {

            // Initialise computation
            register int jpar    = inx[jdev];
            double       g       = (*m_wrk_grad)(jpar);
            double       fa_i    = fa * g;

            // Update gradient.
            (*m_gradient)(jpar) -= fb * g;

            // Loop over rows
            register int* ipar = inx;
            for (register int idev = 0; idev < ndev; ++idev, ++ipar)
                values[idev] = fa_i * (*m_wrk_grad)(*ipar);

            // Add column to matrix
            m_covar->add_col(values, inx, ndev, jpar);

        } // endfor: looped over columns

    } // endfor: iterated over all events

    // Free temporary memory
    if (values != NULL) delete [] values;
    if (inx    != NULL) delete [] inx;

    // Optionally dump gradient and covariance matrix
    #if G_EVAL_DEBUG
    std::cout << *m_gradient << std::endl;
    std::cout << *m_covar << std::endl;
    #endif

    // Timing measurement
    #if G_EVAL_TIMING
    double t_elapse = (double)(clock() - t_start) / (double)CLOCKS_PER_SEC;
    std::cout << "GObservations::optimizer::poisson_unbinned: CPU usage = "
              << t_elapse << " sec" << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Evaluate log-likelihood function for Poisson statistics and
 * binned analysis
 *
 * @param[in] obs Observation.
 * @param[in] pars Optimizer parameters.
 *
 * Evaluate log-likelihood function, gradient and curvature matrix for one
 * observation using binned analysis and Poisson statistics.
 *
 * @todo We still use the pointing information from the event. This has to
 * be removed later ...
 * @todo Need to implement methods that get binsize.
 ***************************************************************************/
void GObservations::optimizer::poisson_binned(const GObservation& obs,
                                              const GOptimizerPars& pars) 
{
    // Timing measurement
    #if G_EVAL_TIMING
    clock_t t_start = clock();
    #endif

    // Get number of parameters
    int npars = pars.npars();

    // Allocate some working arrays
    int*    inx    = new int[npars];
    double* values = new double[npars];

    // Iterate over all bins
    for (int i = 0; i < obs.events()->size(); ++i) {

        // Get pointer to bin
        GEvent* bin = obs.events()->pointer(i);

        // Get number of counts in bin
        double data = bin->counts();

        // Get model and derivative
        double model = obs.model((GModels&)pars, *(bin->dir()),
                                 *(bin->energy()), *(bin->time()),
                                 *(bin->pnt()), m_wrk_grad);

        // Multiply by bin size        
        model *= bin->size();

        // Skip bin if model is too small (avoids -Inf or NaN gradients)
        if (model <= m_minmod)
            continue;

        // Create index array of non-zero derivatives
        int ndev = 0;
        for (int i = 0; i < npars; ++i) {
            if ((*m_wrk_grad)(i) != 0.0 && !std::isinf((*m_wrk_grad)(i))) {
                inx[ndev] = i;
                ndev++;
            }
        }

        // Update Poissonian statistics (excluding factorial term for faster
        // computation)
        m_value -= data * log(model) - model;

        // Skip bin now if there are no non-zero derivatives
        if (ndev < 1)
            continue;

        // Update gradient vector and curvature matrix. To avoid
        // unneccessary computations we distinguish the case where
        // data>0 and data=0. The second case requires much less
        // computation since it does not contribute to the covariance
        // matrix ...
        if (data > 0.0) {

            // Pre computation
            double fb = data / model;
            double fc = (1.0 - fb);
            double fa = fb / model;

            // Loop over columns
            for (int jdev = 0; jdev < ndev; ++jdev) {

                // Initialise computation
                register int jpar    = inx[jdev];
                double       g       = (*m_wrk_grad)(jpar);
                double       fa_i    = fa * g;

                // Update gradient
                (*m_gradient)(jpar) += fc * g;

                // Loop over rows
                register int* ipar = inx;
                for (register int idev = 0; idev < ndev; ++idev, ++ipar)
                    values[idev] = fa_i * (*m_wrk_grad)(*ipar);

                // Add column to matrix
                m_covar->add_col(values, inx, ndev, jpar);
            } // endfor: looped over columns
        } // endif: data was > 0

        // ... handle now data=0
        else {
            register int* ipar = inx;
            for (register int idev = 0; idev < ndev; ++idev, ++ipar)
                (*m_gradient)(*ipar) += (*m_wrk_grad)(*ipar);
        }

    } // endfor: iterated over all events

    // Free temporary memory
    if (values != NULL) delete [] values;
    if (inx    != NULL) delete [] inx;

    // Optionally dump gradient and covariance matrix
    #if G_EVAL_DEBUG
    std::cout << *m_gradient << std::endl;
    std::cout << *m_covar << std::endl;
    #endif

    // Timing measurement
    #if G_EVAL_TIMING
    double t_elapse = (double)(clock() - t_start) / (double)CLOCKS_PER_SEC;
    std::cout << "GObservations::optimizer::poisson_binned: CPU usage = "
              << t_elapse << " sec" << std::endl;
    #endif

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GObservations::optimizer::init_members(void)
{
    // Initialise members
    m_value      = 0.0;
    m_npred      = 0.0;
    m_minmod     = 1.0e-100;
    m_gradient   = NULL;
    m_covar      = NULL;
    m_this       = NULL;
    m_wrk_grad   = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] fct Members to be copied.
 ***************************************************************************/
void GObservations::optimizer::copy_members(const optimizer& fct)
{
    // Copy attributes
    m_value  = fct.m_value;
    m_npred  = fct.m_npred;
    m_minmod = fct.m_minmod;

    // Copy gradient if it exists
    if (fct.m_gradient != NULL)
        m_gradient = new GVector(*fct.m_gradient);

    // Copy covariance matrix if it exists
    if (fct.m_covar != NULL)
        m_covar = new GSparseMatrix(*fct.m_covar);

    // Copy working gradient if it exists
    if (fct.m_wrk_grad != NULL)
        m_wrk_grad = new GVector(*fct.m_wrk_grad);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GObservations::optimizer::free_members(void)
{
    // Free members
    if (m_gradient != NULL) delete m_gradient;
    if (m_covar    != NULL) delete m_covar;
    if (m_wrk_grad != NULL) delete m_wrk_grad;

    // Signal free pointers
    m_gradient = NULL;
    m_covar    = NULL;
    m_wrk_grad = NULL;

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                  Friends                                =
 =                                                                         =
 ==========================================================================*/
