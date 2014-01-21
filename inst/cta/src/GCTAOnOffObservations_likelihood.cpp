/***************************************************************************
 *    GCTAOnOffObservations_likelihood.cpp - Likelihood function class     *
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
 * @file GCTAOnOffObservations_likelihood.cpp
 * @brief Likelihood function class implementation
 * @author Pierrick Martin
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GCTAOnOffObservations.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_EVAL    "GCTAOnOffObservations::likelihood::eval(GOptimizerPars&)"

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
GCTAOnOffObservations::likelihood::likelihood(void) : GOptimizerFunction()
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
 * Constructs optimizer from GCTAOnOffObservations container. The method copies the
 * pointer to the observation container in the m_this member, making the
 * observation container accessible to the optimizer class.
 ***************************************************************************/
GCTAOnOffObservations::likelihood::likelihood(GCTAOnOffObservations *obs) : GOptimizerFunction()
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
GCTAOnOffObservations::likelihood::likelihood(const likelihood& fct) :
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
GCTAOnOffObservations::likelihood::~likelihood(void)
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
GCTAOnOffObservations::likelihood& GCTAOnOffObservations::likelihood::operator=(const likelihood& fct)
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
 * optimization. It handles ON-OFF analysis. The likelihood function and
 * the derivatives are computed in GCTAOnOffObservation.poisson_onoff().
 ***************************************************************************/
void GCTAOnOffObservations::likelihood::eval(const GOptimizerPars& pars) 
{
    // Timing measurement
    #if defined(G_EVAL_TIMING)
    clock_t t_start = clock();
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
        if (m_gradient  != NULL) delete m_gradient;
        if (m_curvature != NULL) delete m_curvature;

        // Initialise value, gradient vector and curvature matrix
        m_value     = 0.0;
        m_npred     = 0.0;
        m_gradient  = new GVector(npars);
        m_curvature = new GMatrixSparse(npars,npars);
		
		// Loop over observations in container
		double v=0.0;
		for (int i = 0; i < m_this->size(); ++i) {
			
			// Compute likelihood
			v += m_this->m_obs[i]->likelihood_poisson_onoff(pars,
											           m_curvature,
											           m_gradient,
											           m_value,
											           m_npred);
			
		} // endfor: looped over observations	

    } while(0); // endwhile: main loop

    // Copy over the parameter gradients for all parameters that are
    // free (so that we can access the gradients from outside)
    for (int ipar = 0; ipar < pars.size(); ++ipar) {
        if (pars[ipar]->is_free()) {
            GOptimizerPar* par = const_cast<GOptimizerPar*>(pars[ipar]);
            par->factor_gradient((*m_gradient)[ipar]);
        }
    }

    // Optionally dump gradient and curvature matrix
    #if defined(G_EVAL_DEBUG)
    std::cout << *m_gradient << std::endl;
    for (int i = 0; i < pars.size(); ++i) {
        for (int j = 0; j < pars.size(); ++j) {
            std::cout << (*m_curvature)(i,j) << " ";
        }
        std::cout << std::endl;
    }
    #endif

    // Timing measurement
    #if defined(G_EVAL_TIMING)
    double t_elapse = (double)(clock() - t_start) / (double)CLOCKS_PER_SEC;
    std::cout << "GCTAOnOffObservations::optimizer::eval: CPU usage = "
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
void GCTAOnOffObservations::likelihood::init_members(void)
{
    // Initialise members
    m_value     = 0.0;
    m_npred     = 0.0;
    m_this      = NULL;
    m_gradient  = NULL;
    m_curvature = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] fct Optimizer.
 ***************************************************************************/
void GCTAOnOffObservations::likelihood::copy_members(const likelihood& fct)
{
    // Copy attributes
    m_value  = fct.m_value;
    m_npred  = fct.m_npred;
    m_this   = fct.m_this;

    // Clone gradient if it exists
    if (fct.m_gradient != NULL) m_gradient = new GVector(*fct.m_gradient);

    // Clone curvature matrix if it exists
    if (fct.m_curvature != NULL) m_curvature = new GMatrixSparse(*fct.m_curvature);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAOnOffObservations::likelihood::free_members(void)
{
    // Free members
    if (m_gradient  != NULL) delete m_gradient;
    if (m_curvature != NULL) delete m_curvature;

    // Signal free pointers
    m_gradient  = NULL;
    m_curvature = NULL;

    // Return
    return;
}
