/***************************************************************************
 *  GObservations_optimizer.cpp  -  Optimizer class of observations class  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2012 by Juergen Knoedlseder                         *
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
 * @file GObservations_optimizer.cpp
 * @brief Model parameter optimization class implementation
 * @author J. Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GObservations.hpp"
#include "GTools.hpp"
#include "GEvent.hpp"
#include "GEventList.hpp"
#include "GEventCube.hpp"
#include "GEventBin.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_EVAL              "GObservations::optimizer::eval(GOptimizerPars&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
#define G_EVAL_TIMING   0 //!< Perform optimizer timing (0=no, 1=yes)
#define G_EVAL_DEBUG    0 //!< Perform optimizer debugging (0=no, 1=yes)
#define G_OPT_DEBUG     0 //!< Perform optimizer debugging (0=no, 1=yes)

/* __ Prototypes _________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
***************************************************************************/
GObservations::optimizer::optimizer(void) : GOptimizerFunction()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Construct optimizer from GObservations observation container
 *
 * @param[in] obs Observations container.
 ***************************************************************************/
GObservations::optimizer::optimizer(GObservations *obs) : GOptimizerFunction()
{
    // Initialise members
    init_members();

    // Set object
    m_this = obs;

    // Return
    return;
}



/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] fct Optimizer.
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
GObservations::optimizer::~optimizer(void)
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
 * @param[in] fct Optimizer.
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
 * @exception GException::invalid_statistics
 *            Invalid optimization statistics encountered.
 *
 * This method evaluates the -(log-likelihood) function for parameter
 * optimisation. It handles both binned and unbinned data and supportes
 * Poisson and Gaussian statistics. 
 * Note that different statistics and different analysis methods
 * (binned/unbinned) may be combined.
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
        if (npars < 1) {
            continue;
        }

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
        
        //vector to save working variables of each thread.
        std::vector<GVector* > vect_cpy_mgrad;
        std::vector<GSparseMatrix* > vect_cpy_mcovar;
        std::vector<double* > vect_cpy_mvalue;
        std::vector<double* > vect_cpy_mnpred;
        
        // Here OpenMP will paralellize the execution. The following code will be executed by the differents threads.
        // In order to avoid protecting attributes ( m_value,m_npred, m_gradient and m_covar), each thread works with its own working variables (cpy_*)
        // When a thread starts, we add working variables in a vector (vect_cpy_*). When computation is finished we just add all elements contain in the vector to the attributes value.
        #pragma omp parallel
        {
            // Copy variables for multi-threading
            GModels cpy_model((GModels&)pars);
            
            GVector cpy_wrk_grad(npars);
            GVector* cpy_mgrad = new GVector(npars);
            GSparseMatrix* cpy_mcovar = new GSparseMatrix(npars,npars);
            cpy_mcovar->stack_init(npars,10000);
        
            double* cpy_mnpred = new double;
            double* cpy_mvalue = new double;
            
            //Critical zone
            #pragma omp critical
            {
                vect_cpy_mgrad.push_back(cpy_mgrad);
                vect_cpy_mcovar.push_back(cpy_mcovar); 
                vect_cpy_mvalue.push_back(cpy_mvalue);
                vect_cpy_mnpred.push_back(cpy_mnpred);
            }

            // The omp for directive will deal the iterations on the differents threads.
            #pragma omp for
            // Loop over all observations
            for (int i = 0; i < m_this->size(); ++i) {
    
                // Extract statistics for this observation
                std::string statistics = m_this->m_obs[i]->statistics();
    
                // Unbinned analysis
                if (dynamic_cast<const GEventList*>(m_this->m_obs[i]->events()) != NULL) {
    
                    // Poisson statistics
                    if (toupper(statistics) == "POISSON") {
    
                        // Determine Npred value and gradient for this observation
                        double npred = m_this->m_obs[i]->npred(cpy_model, &cpy_wrk_grad);
    
                        // Update the Npred value, gradient.
                        *cpy_mnpred += npred;
                        *cpy_mgrad+=cpy_wrk_grad;

                        #if G_EVAL_DEBUG
                        #pragma omp critial single
                        {
                            std::cout << "Unbinned Poisson (" << i << "):";
                            std::cout << " Npred=" << npred;
                            std::cout << " Grad="<< *copy_wrk_grad << std::endl;
                            std::cout << "Sum:";
                            std::cout << " Npred=" << *cpy_mnpred;
                            std::cout << " Grad="<< *cpy_mgrad<< std::endl;
                        }
                        #endif
    
                        // Update the log-likelihood
                        poisson_unbinned(*(m_this->m_obs[i]), cpy_model,*cpy_mcovar,*cpy_mgrad,*cpy_mvalue,cpy_wrk_grad);
    
                        // Add the Npred value to the log-likelihood
                        *cpy_mvalue += npred;
    
                    } // endif: Poisson statistics
    
                    // ... otherwise throw an exception
                    else
                        throw GException::invalid_statistics(G_EVAL, statistics,
                            "Unbinned optimization requires Poisson statistics.");
    
                } // endif: unbinned analysis
    
                // ... or binned analysis
                else {
    
                    // Poisson statistics
                    if (toupper(statistics) == "POISSON") {
                        #if G_EVAL_DEBUG
                        std::cout << "Binned Poisson" << std::endl;
                        #endif
                        poisson_binned(*(m_this->m_obs[i]), cpy_model,*cpy_mcovar,*cpy_mgrad,*cpy_mvalue,*cpy_mnpred,cpy_wrk_grad);
                    }
    
                    // ... or Gaussian statistics
                    else if (toupper(statistics) == "GAUSSIAN") {
                        #if G_EVAL_DEBUG
                        std::cout << "Binned Gaussian" << std::endl;
                        #endif
                        gaussian_binned(*(m_this->m_obs[i]), pars);
                    }
    
                    // ... or unsupported
                    else
                        throw GException::invalid_statistics(G_EVAL, statistics,
                            "Binned optimization requires Poisson or Gaussian statistics.");
    
                } // endelse: binned analysis
    
            } // endfor: looped over observations
            
        } // end pragma omp parallel
        
        // Now the computation is finished, update attributes.
        // For each omp section, a thread will be created.
        #pragma omp sections
        {
            #pragma omp section
            {
                for(int i =0;i<vect_cpy_mcovar.size();i++){
                    *m_covar += *(vect_cpy_mcovar.at(i));
                    delete vect_cpy_mcovar.at(i);
                }
            }

            #pragma omp section
            {
                for(int i =0;i<vect_cpy_mgrad.size();i++){
                    *m_gradient += *(vect_cpy_mgrad.at(i));
                    delete vect_cpy_mgrad.at(i);
                }
            }

            #pragma omp section
            {
                for(int i =0;i<vect_cpy_mnpred.size();i++){
                    m_npred += *(vect_cpy_mnpred.at(i));
                    delete vect_cpy_mnpred.at(i);
                }
            }

            #pragma omp section
            {
                for(int i =0;i<vect_cpy_mvalue.size();i++){
                    m_value += *(vect_cpy_mvalue.at(i));
                    delete vect_cpy_mvalue.at(i);
                }
            }
    }

        // Release stack
        m_covar->stack_destroy();

    } while(0); // endwhile: main loop
    
    // Copy over the parameter gradients for all parameters that are
    // free (so that we can access the gradients from outside)
    for (int ipar = 0; ipar < pars.npars(); ++ipar) {
        if (pars.par(ipar).isfree()) {
            GModelPar* par = (GModelPar*)&pars.par(ipar);
            par->gradient((*m_gradient)[ipar]);
        }
    }

    // Optionally dump gradient and covariance matrix
    #if G_EVAL_DEBUG
    std::cout << *m_gradient << std::endl;
    for (int i = 0; i < pars.npars(); ++i) {
        for (int j = 0; j < pars.npars(); ++j) {
            std::cout << (*m_covar)(i,j) << " ";
        }
        std::cout << std::endl;
    }
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
 *        unbinned analysis
 *
 * @param[in] obs Observation.
 * @param[in] pars Optimizer parameters.
 *
 * This method evaluates the -(log-likelihood) function for parameter
 * optimisation using unbinned analysis and Poisson statistics.
 * The -(log-likelihood) function is given by
 * \f$L=-\sum_i \log e_i\f$
 * where the sum is taken over all events. This method also computes the
 * parameter gradients
 * \f$\delta L/dp\f$
 * and the curvature matrix
 * \f$\delta^2 L/dp_1 dp_2\f$.
 ***************************************************************************/
void GObservations::optimizer::poisson_unbinned(const GObservation& obs,
                                                const GOptimizerPars& pars)
{
    poisson_unbinned(obs,pars,*m_covar,*m_gradient,m_value,*m_wrk_grad); 
}

void GObservations::optimizer::poisson_unbinned(const GObservation& obs,
                                                const GOptimizerPars& pars,
                                                GSparseMatrix& covar,
                                                GVector& mgrad,
                                                double& value,
                                                GVector& gradient)
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

        // Get event pointer
        const GEvent* event = (*obs.events())[i];

        // Get model and derivative
        double model = obs.model((GModels&)pars, *event, &gradient);

        // Skip bin if model is too small (avoids -Inf or NaN gradients)
        if (model <= m_minmod) {
            continue;
        }

        // Create index array of non-zero derivatives and initialise working
        // array
        int ndev = 0;
        for (int i = 0; i < npars; ++i) {
            values[i] = 0.0;
            if (gradient[i] != 0.0 && !isinfinite(gradient[i])) {
                inx[ndev] = i;
                ndev++;
            }
        }

        // Update Poissonian statistics (excluding factorial term for faster
        // computation)
        value -= log(model);

        // Skip bin now if there are no non-zero derivatives
        if (ndev < 1) {
            continue;
        }

        // Update gradient vector and curvature matrix.
        double fb = 1.0 / model;
        double fa = fb / model;
        for (int jdev = 0; jdev < ndev; ++jdev) {

            // Initialise computation
            register int jpar    = inx[jdev];
            double       g       = gradient[jpar];
            double       fa_i    = fa * g;

            // Update gradient.
            mgrad[jpar] -= fb * g;

            // Loop over rows
            register int* ipar = inx;

            for (register int idev = 0; idev < ndev; ++idev, ++ipar) {
                values[idev] = fa_i * gradient[*ipar];
            }

            // Add column to matrix
            covar.add_col(values, inx, ndev, jpar);

        } // endfor: looped over columns

    } // endfor: iterated over all events

    // Free temporary memory
    if (values != NULL) delete [] values;
    if (inx    != NULL) delete [] inx;

    // Optionally dump gradient and covariance matrix
    #if G_EVAL_DEBUG
    std::cout << mgrad << std::endl;
    std::cout << covar << std::endl;
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
 *        binned analysis
 *
 * @param[in] obs Observation.
 * @param[in] pars Optimizer parameters.
 *
 * This method evaluates the -(log-likelihood) function for parameter
 * optimisation using binned analysis and Poisson statistics.
 * The -(log-likelihood) function is given by
 * \f$L=-\sum_i n_i \log e_i - e_i\f$
 * where the sum is taken over all data space bins, \f$n_i\f$ is the
 * observed number of counts and \f$e_i\f$ is the model.
 * This method also computes the parameter gradients
 * \f$\delta L/dp\f$
 * and the curvature matrix
 * \f$\delta^2 L/dp_1 dp_2\f$
 * and also updates the total number of predicted events m_npred.
 ***************************************************************************/
void GObservations::optimizer::poisson_binned(const GObservation& obs,
                                              const GOptimizerPars& pars) 
{
    poisson_binned(obs,pars,*m_covar,*m_gradient,m_value,m_npred,*m_wrk_grad); 
}

void GObservations::optimizer::poisson_binned(const GObservation& obs,
                                                const GOptimizerPars& pars,
                                                GSparseMatrix& covar,
                                                GVector& mgrad,
                                                double& value,
                                                double& npred,
                                                GVector& gradient)
{
    // Timing measurement
    #if G_EVAL_TIMING
    clock_t t_start = clock();
    #endif

    // Initialise statistics
    #if G_OPT_DEBUG
    int    n_bins        = 0;
    int    n_used        = 0;
    int    n_small_model = 0;
    int    n_zero_data   = 0;
    double sum_data      = 0.0;
    double sum_model     = 0.0;
    double init_value    = value;
    #endif

    // Get number of parameters
    int npars = pars.npars();

    // Allocate some working arrays
    int*    inx    = new int[npars];
    double* values = new double[npars];

    // Iterate over all bins
    for (int i = 0; i < obs.events()->size(); ++i) {

        // Update number of bins
        #if G_OPT_DEBUG
        n_bins++;
        #endif

        // Get event pointer
        const GEventBin* bin = (*((GEventCube*)obs.events()))[i];

        // Get number of counts in bin
        double data = bin->counts();

        // Get model and derivative
        double model = obs.model((GModels&)pars, *bin, &gradient);

        // Multiply model by bin size
        model *= bin->size();

        // Skip bin if model is too small (avoids -Inf or NaN gradients)
        if (model <= m_minmod) {
            #if G_OPT_DEBUG
            n_small_model++;
            #endif
            continue;
        }

        // Update statistics
        #if G_OPT_DEBUG
        n_used++;
        sum_data  += data;
        sum_model += model;
        #endif

        // Update Npred
        npred += model;

        // Multiply gradient by bin size
        gradient *= bin->size();

        // Create index array of non-zero derivatives and initialise working
        // array
        int ndev = 0;
        for (int i = 0; i < npars; ++i) {
            values[i] = 0.0;
            if (gradient[i] != 0.0 && !isinfinite(gradient[i])) {
                inx[ndev] = i;
                ndev++;
            }
        }

        // Update gradient vector and curvature matrix. To avoid
        // unneccessary computations we distinguish the case where
        // data>0 and data=0. The second case requires much less
        // computation since it does not contribute to the covariance
        // matrix ...
        if (data > 0.0) {

            // Update Poissonian statistics (excluding factorial term for
            // faster computation)
            value -= data * log(model) - model;

            // Skip bin now if there are no non-zero derivatives
            if (ndev < 1) {
                continue;
            }

            // Pre computation
            double fb = data / model;
            double fc = (1.0 - fb);
            double fa = fb / model;

            // Loop over columns
            for (int jdev = 0; jdev < ndev; ++jdev) {

                // Initialise computation
                register int jpar    = inx[jdev];
                double       g       = gradient[jpar];
                double       fa_i    = fa * g;

                // Update gradient
                (*m_gradient)[jpar] += fc * g;

                // Loop over rows
                register int* ipar = inx;
                for (register int idev = 0; idev < ndev; ++idev, ++ipar) {
                    values[idev] = fa_i * gradient[*ipar];
                }

                // Add column to matrix
                covar.add_col(values, inx, ndev, jpar);

            } // endfor: looped over columns

        } // endif: data was > 0

        // ... handle now data=0
        else {

            // Update statistics
            #if G_OPT_DEBUG
            n_zero_data++;
            #endif
            
            // Update Poissonian statistics (excluding factorial term for
            // faster computation)
            value += model;

            // Skip bin now if there are no non-zero derivatives
            if (ndev < 1) {
                continue;
            }

            // Update gradient
            register int* ipar = inx;
            for (register int idev = 0; idev < ndev; ++idev, ++ipar) {
                mgrad[*ipar] += gradient[*ipar];
            }

        } // endif: data was 0

    } // endfor: iterated over all events

    // Free temporary memory
    if (values != NULL) delete [] values;
    if (inx    != NULL) delete [] inx;

    // Dump statistics
    #if G_OPT_DEBUG
    std::cout << "Number of bins: " << n_bins << std::endl;
    std::cout << "Number of bins used for computation: " << n_used << std::endl;
    std::cout << "Number of bins excluded due to small model: " << n_small_model << std::endl;
    std::cout << "Number of bins with zero data: " << n_zero_data << std::endl;
    std::cout << "Sum of data: " << sum_data << std::endl;
    std::cout << "Sum of model: " << sum_model << std::endl;
    std::cout << "Initial statistics: " << init_value << std::endl;
    std::cout << "Statistics: " << m_value-init_value << std::endl;
    #endif

    // Optionally dump gradient and covariance matrix
    #if G_EVAL_DEBUG
    std::cout << mgrad << std::endl;
    std::cout << covar << std::endl;
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


/***********************************************************************//**
 * @brief Evaluate log-likelihood function for Gaussian statistics and
 *        binned analysis
 *
 * @param[in] obs Observation.
 * @param[in] pars Optimizer parameters.
 *
 * This method evaluates the -(log-likelihood) function for parameter
 * optimisation using binned analysis and Poisson statistics.
 * The -(log-likelihood) function is given by
 * \f$L = 1/2 \sum_i (n_i - e_i)^2 \sigma_i^{-2}\f$
 * where the sum is taken over all data space bins, \f$n_i\f$ is the
 * observed number of counts, \f$e_i\f$ is the model and \f$\sigma_i\f$
 * is the statistical uncertainty.
 * This method also computes the parameter gradients
 * \f$\delta L/dp\f$
 * and the curvature matrix
 * \f$\delta^2 L/dp_1 dp_2\f$
 * and also updates the total number of predicted events m_npred.
 ***************************************************************************/
void GObservations::optimizer::gaussian_binned(const GObservation& obs,
                                               const GOptimizerPars& pars) 
{
    gaussian_binned(obs,pars,*m_covar,*m_gradient,m_value,m_npred,*m_wrk_grad); 
}
void GObservations::optimizer::gaussian_binned(const GObservation& obs,
                                              const GOptimizerPars& pars,
                                              GSparseMatrix& covar,
                                              GVector& mgrad,
                                              double& value,
                                              double& npred,
                                              GVector& gradient)
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

        // Get event pointer
        const GEventBin* bin = (*((GEventCube*)obs.events()))[i];

        // Get number of counts in bin
        double data = bin->counts();

        // Get statistical uncertainty
        double sigma = bin->error();

        // Skip bin if statistical uncertainty is too small
        if (sigma <= m_minerr) {
            continue;
        }

        // Get model and derivative
        double model = obs.model((GModels&)pars, *bin, &gradient);

        // Multiply model by bin size
        model *= bin->size();

        // Skip bin if model is too small (avoids -Inf or NaN gradients)
        if (model <= m_minmod) {
            continue;
        }

        // Update Npred
        npred += model;

        // Multiply gradient by bin size
        gradient *= bin->size();

        // Create index array of non-zero derivatives and initialise working
        // array
        int ndev = 0;
        for (int i = 0; i < npars; ++i) {
            values[i] = 0.0;
            if (gradient[i] != 0.0 && !isinfinite(gradient[i])) {
                inx[ndev] = i;
                ndev++;
            }
        }

        // Set weight
        double weight = 1.0 / (sigma * sigma);

        // Update Gaussian statistics
        double fa = data - model;
        value  += 0.5 * (fa * fa * weight);
            
        // Skip bin now if there are no non-zero derivatives
        if (ndev < 1) {
            continue;
        }

        // Loop over columns
        for (int jdev = 0; jdev < ndev; ++jdev) {

            // Initialise computation
            register int jpar = inx[jdev];
            double       fa_i = gradient[jpar] * weight;

            // Update gradient
            gradient[jpar] -= fa * fa_i;

            // Loop over rows
            register int* ipar = inx;
            for (register int idev = 0; idev < ndev; ++idev, ++ipar) {
                values[idev] = fa_i * gradient[*ipar];
            }

            // Add column to matrix
            covar.add_col(values, inx, ndev, jpar);

        } // endfor: looped over columns

    } // endfor: iterated over all events

    // Free temporary memory
    if (values != NULL) delete [] values;
    if (inx    != NULL) delete [] inx;

    // Optionally dump gradient and covariance matrix
    #if G_EVAL_DEBUG
    std::cout << gradient << std::endl;
    std::cout << covar << std::endl;
    #endif

    // Timing measurement
    #if G_EVAL_TIMING
    double t_elapse = (double)(clock() - t_start) / (double)CLOCKS_PER_SEC;
    std::cout << "GObservations::optimizer::gaussian_binned: CPU usage = "
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
    m_minerr     = 1.0e-100;
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
 * @param[in] fct Optimizer.
 ***************************************************************************/
void GObservations::optimizer::copy_members(const optimizer& fct)
{
    // Copy attributes
    m_value  = fct.m_value;
    m_npred  = fct.m_npred;
    m_minmod = fct.m_minmod;
    m_minerr = fct.m_minerr;

    // Clone gradient if it exists
    if (fct.m_gradient != NULL) m_gradient = new GVector(*fct.m_gradient);

    // Clone covariance matrix if it exists
    if (fct.m_covar != NULL) m_covar = new GSparseMatrix(*fct.m_covar);

    // Clone working gradient if it exists
    if (fct.m_wrk_grad != NULL) m_wrk_grad = new GVector(*fct.m_wrk_grad);

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
