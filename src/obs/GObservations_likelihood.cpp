/***************************************************************************
 *         GObservations_likelihood.cpp - Likelihood function class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2013 by Juergen Knoedlseder                         *
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
 * @file GObservations_likelihood.cpp
 * @brief Likelihood function class implementation
 * @author Juergen Knoedlseder
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

/* __ OpenMP section _____________________________________________________ */
#ifdef _OPENMP
#include <omp.h>
#endif

/* __ Method name definitions ____________________________________________ */
#define G_EVAL             "GObservations::likelihood::eval(GOptimizerPars&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_EVAL_TIMING     //!< Perform optimizer timing (0=no, 1=yes)
//#define G_EVAL_DEBUG      //!< Perform optimizer debugging (0=no, 1=yes)

/* __ Prototypes _________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GObservations::likelihood::likelihood(void) : GOptimizerFunction()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Observations constructor
 *
 * @param[in] obs Observations container pointer.
 *
 * Constructs optimizer from GObservations container. The method copies the
 * pointer to the observation container in the m_this member, making the
 * observation container accessible to the optimizer class.
 ***************************************************************************/
GObservations::likelihood::likelihood(GObservations *obs) : GOptimizerFunction()
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
 * @param[in] fct Optimizer function.
 ***************************************************************************/
GObservations::likelihood::likelihood(const likelihood& fct) :
                           GOptimizerFunction(fct)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(fct);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GObservations::likelihood::~likelihood(void)
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
 * @param[in] fct Likelihood function.
 * @return Likelihood function.
 ***************************************************************************/
GObservations::likelihood& GObservations::likelihood::operator=(const likelihood& fct)
{
    // Execute only if object is not identical
    if (this != &fct) {

        // Copy base class members
        this->GOptimizerFunction::operator=(fct);

        // Free members
        free_members();

        // Initialise private members
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
 * optimization. It handles both binned and unbinned data and supportes
 * Poisson and Gaussian statistics. 
 * Note that different statistics and different analysis methods
 * (binned/unbinned) may be combined.
 ***************************************************************************/
void GObservations::likelihood::eval(const GOptimizerPars& pars) 
{
    // Timing measurement
    #if defined(G_EVAL_TIMING)
    #ifdef _OPENMP
    double t_start = omp_get_wtime();
    #else
    clock_t t_start = clock();
    #endif
    #endif

    // Single loop for common exit point
    do {

        // Get number of parameters
        int npars = pars.size();

        // Fall through if we have no free parameters
        if (npars < 1) {
            continue;
        }

        // Free old memory
        if (m_gradient != NULL) delete m_gradient;
        if (m_covar    != NULL) delete m_covar;

        // Initialise value, gradient vector and curvature matrix
        m_value    = 0.0;
        m_npred    = 0.0;
        m_gradient = new GVector(npars);
        m_covar    = new GMatrixSparse(npars,npars);

        // Set stack size and number of entries
        int stack_size  = (2*npars > 100000) ? 2*npars : 100000;
        int max_entries =  2*npars;
        m_covar->stack_init(stack_size, max_entries);

        // Allocate vectors to save working variables of each thread
        std::vector<GVector*>       vect_cpy_grad;
        std::vector<GMatrixSparse*> vect_cpy_covar;
        std::vector<double*>        vect_cpy_value;
        std::vector<double*>        vect_cpy_npred;

        // Here OpenMP will paralellize the execution. The following code will
        // be executed by the differents threads. In order to avoid protecting
        // attributes ( m_value,m_npred, m_gradient and m_covar), each thread
        // works with its own working variables (cpy_*). When a thread starts,
        // we add working variables in a vector (vect_cpy_*). When computation
        // is finished we just add all elements contain in the vector to the
        // attributes value.
        #pragma omp parallel
        {
            // Allocate and initialize variable copies for multi-threading
            GModels        cpy_model(m_this->models());
            GVector*       cpy_gradient = new GVector(npars);
            GMatrixSparse* cpy_covar    = new GMatrixSparse(npars,npars);
            double*        cpy_npred    = new double(0.0);
            double*        cpy_value    = new double(0.0);

            // Set stack size and number of entries
            cpy_covar->stack_init(stack_size, max_entries);

            // Push variable copies into vector. This is a critical zone to
            // avoid multiple thread pushing simultaneously.
            #pragma omp critical
            {
                vect_cpy_grad.push_back(cpy_gradient);
                vect_cpy_covar.push_back(cpy_covar); 
                vect_cpy_value.push_back(cpy_value);
                vect_cpy_npred.push_back(cpy_npred);
            }

            // Loop over all observations. The omp for directive will deal
            // with the iterations on the differents threads.
            #pragma omp for
            for (int i = 0; i < m_this->size(); ++i) {

                // Compute likelihood
                *cpy_value += m_this->m_obs[i]->likelihood(cpy_model,
                                                           cpy_gradient,
                                                           cpy_covar,
                                                           cpy_npred);

            } // endfor: looped over observations

            // Release stack
            cpy_covar->stack_destroy();

        } // end pragma omp parallel

        // Now the computation is finished, update attributes.
        // For each omp section, a thread will be created.
        #pragma omp sections
        {
            #pragma omp section
            {
                for (int i = 0; i < vect_cpy_covar.size() ; ++i) {
                    *m_covar += *(vect_cpy_covar.at(i));
                    delete vect_cpy_covar.at(i);
                }
            }

            #pragma omp section
            {
                for (int i = 0; i < vect_cpy_grad.size(); ++i){
                    *m_gradient += *(vect_cpy_grad.at(i));
                    delete vect_cpy_grad.at(i);
                }
            }

            #pragma omp section
            {
                for(int i = 0; i < vect_cpy_npred.size(); ++i){
                    m_npred += *(vect_cpy_npred.at(i));
                    delete vect_cpy_npred.at(i);
                }
            }

            #pragma omp section
            {
                for (int i = 0; i < vect_cpy_value.size(); ++i){
                    m_value += *(vect_cpy_value.at(i));
                    delete vect_cpy_value.at(i);
                }
            }
        } // end of pragma omp sections

        // Release stack
        m_covar->stack_destroy();

    } while(0); // endwhile: main loop

    // Copy over the parameter gradients for all parameters that are
    // free (so that we can access the gradients from outside)
    for (int ipar = 0; ipar < pars.size(); ++ipar) {
        if (pars[ipar]->isfree()) {
            GOptimizerPar* par = const_cast<GOptimizerPar*>(pars[ipar]);
            par->factor_gradient((*m_gradient)[ipar]);
        }
    }

    // Optionally dump gradient and covariance matrix
    #if defined(G_EVAL_DEBUG)
    std::cout << *m_gradient << std::endl;
    for (int i = 0; i < pars.size(); ++i) {
        for (int j = 0; j < pars.size(); ++j) {
            std::cout << (*m_covar)(i,j) << " ";
        }
        std::cout << std::endl;
    }
    #endif

    // Timing measurement
    #if defined(G_EVAL_TIMING)
    #ifdef _OPENMP
    double t_elapse = omp_get_wtime()-t_start;
    #else
    double t_elapse = (double)(clock() - t_start) / (double)CLOCKS_PER_SEC;
    #endif
    std::cout << "GObservations::optimizer::eval: CPU usage = "
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
void GObservations::likelihood::init_members(void)
{
    // Initialise members
    m_value     = 0.0;
    m_npred     = 0.0;
    m_this      = NULL;
    m_gradient  = NULL;
    m_covar     = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] fct Optimizer.
 ***************************************************************************/
void GObservations::likelihood::copy_members(const likelihood& fct)
{
    // Copy attributes
    m_value  = fct.m_value;
    m_npred  = fct.m_npred;
    m_this   = fct.m_this;

    // Clone gradient if it exists
    if (fct.m_gradient != NULL) m_gradient = new GVector(*fct.m_gradient);

    // Clone covariance matrix if it exists
    if (fct.m_covar != NULL) m_covar = new GMatrixSparse(*fct.m_covar);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GObservations::likelihood::free_members(void)
{
    // Free members
    if (m_gradient != NULL) delete m_gradient;
    if (m_covar    != NULL) delete m_covar;

    // Signal free pointers
    m_gradient = NULL;
    m_covar    = NULL;

    // Return
    return;
}
