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
		for (int i = 0; i < m_this->size(); ++i) {
			
			// Compute likelihood
			m_value += m_this->m_obs[i]->likelihood_poisson_onoff(m_this->models(),
											                m_curvature,
															m_gradient,
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


/***********************************************************************//**
 * @brief Compute Hessian matrix
 *
 * @param[in] pars Optimizer parameters.
 *
 * @return Hessian matrix.
 ***************************************************************************/
GMatrixSparse GCTAOnOffObservations::likelihood::hessian(const GOptimizerPars& pars)
{
    // Set strategy constants (low)
    //const int    ncyles             = 3;
    //const double step_tolerance     = 0.5;
    //const double gradient_tolerance = 0.1;

    // Set strategy constants (medium)
    const int    ncyles             = 5;
    const double step_tolerance     = 0.3;
    const double gradient_tolerance = 0.05;

    // Set strategy constants (high)
    //const int    ncyles             = 7;
    //const double step_tolerance     = 0.1;
    //const double gradient_tolerance = 0.02;


    // Create working copy of parameters
    GOptimizerPars wrk_pars = pars;

    // Get number of parameters
    int npars = wrk_pars.size();

    // Allocate Hessian matrix
    GMatrixSparse hessian(npars, npars);

    // Find out machine precision
    double eps = 0.1;
    while (1.0+eps != 1.0) {
        eps *= 0.5;
    }
    double eps2 = 2.0 * std::sqrt(eps);

    // Function value
    eval(wrk_pars);
    double f = value();

    // Compute aimsag
    double aimsag = std::sqrt(eps2) * std::abs(f);

    // Diagonal elements
    std::vector<double> g2(npars, 0.0);
    std::vector<double> grd(npars, 0.0);
    std::vector<double> dir(npars, 0.0);
    std::vector<double> yy(npars, 0.0);

    // Loop over parameters
    for (int i = 0; i < npars; ++i) {
    
        // Get parameter
        GOptimizerPar* par = wrk_pars[i];

        // Interrupt if parameter is fixed
        if (par->is_fixed()) {
            hessian(i,i) = 0.0;
            continue;
        }

        // Setup step size
        double xtf  = par->factor_value();
        double dmin = 8.0 * eps2 * std::abs(xtf);
        double d    = 0.000001;
        if (d < dmin) {
            d = dmin;
        }

        // Loop over cycles
        for (int icyc = 0; icyc < ncyles; ++icyc) {
        //for (int icyc = 0; icyc < 1; ++icyc) {

            // Initialise
            double sag = 0.0;
            double fs1 = 0.0; // right-hand side
            double fs2 = 0.0; // left-hand side

            // Compute gradient
            for (int multpy = 0; multpy < 5; ++multpy) {
            //for (int multpy = 0; multpy < 1; ++multpy) {

                // Compute right-hand side
                par->factor_value(xtf + d);
                eval(wrk_pars);
                fs1  = value();
                
                // Compute left-hand side
                par->factor_value(xtf - d);
                eval(wrk_pars);
                fs2  = value();
                
                // Recover current value
                par->factor_value(xtf);
                
                // Compute sag
                sag = 0.5 * (fs1 + fs2 - 2.0*f);

                // Break if sag is okay
                if (std::abs(sag) > eps2 || sag == 0.0) {
                    break;
                }

                // ... otherwise increase step size
                d *= 10.0;
                
            } // endfor

            // Save old step size and second derivative
            double dlast  = d;
            double g2bfor = g2[i];

            // Compute parameter derivatives and store step size and
            // function value
            g2[i]  = 2.0 * sag/(d*d);
            grd[i] = (fs1-fs2)/(2.*d);
            dir[i] = d;
            yy[i]  = fs1;

            // Compute a new step size based on the aimed sag
            if (sag != 0.0) {
                d = std::sqrt(2.0*aimsag/std::abs(g2[i]));
            }
            if (d < dmin) {
                d = dmin;
            }
            /*
            else if (par->factor_value()+d > par->factor_max()) {
                d = dmin;
            }
            else if (par->factor_value()-d > par->factor_min()) {
                d = dmin;
            }
            */
            
            // Check if converged
            if (std::abs((d-dlast)/d) < step_tolerance) {
                break;
            }
            if (std::abs((g2[i]-g2bfor)/g2[i]) < gradient_tolerance) {
                break;
            }
            d = std::min(d, 10.*dlast);
            d = std::max(d, 0.1*dlast);

        } // endfor: cycles

        // Set diagonal element
        hessian(i,i) = g2[i];

    } // endfor: looped over all parameters

    // Debug dump
    #if defined(G_HESSIAN)
    std::cout << "GCTAOnOffObservations::likelihood::hessian: ";
    std::cout << "deltas and gradients:" << std::endl;
    for (int i = 0; i < npars; ++i) {
        std::cout << dir[i] << " ";
        std::cout << grd[i] << std::endl;
    }
    #endif

    // Compute off-diagonal elements
    for (int i = 0; i < npars; ++i) {
    
        // Get parameter 1
        GOptimizerPar* par1 = wrk_pars[i];
        double         x1   = par1->factor_value();
        
        // Increment parameter 1
        par1->factor_value(x1 + dir[i]);

        // Loop over columns
        for (int j = i+1; j < npars; ++j) {

            // Get parameter 2
            GOptimizerPar* par2 = wrk_pars[j];
            double         x2   = par2->factor_value();

            // Interrupt if parameter is fixed
            if (par1->is_fixed() || par2->is_fixed()) {
                hessian(i,j) = 0.0;
                hessian(j,i) = 0.0;
                continue;
            }

            // Increment parameter 2
            par2->factor_value(x2 + dir[j]);

            // Evaluate Hessian element
            eval(wrk_pars);
            double fs1     = value();
            double element = (fs1 + f - yy[i] - yy[j])/(dir[i]*dir[j]);
            
            // Store Hessian element
            hessian(i,j) = element;
            hessian(j,i) = element;

            // Restore parameter 2
            par2->factor_value(x2);
            
        } // endfor: looped over columns
        
        // Restore parameter 1
        par1->factor_value(x1);

    } // endfor: looped over parameters

    // Debug dump
    #if defined(G_HESSIAN)
    std::cout << "GCTAOnOffObservations::likelihood::hessian: " << std::endl;
    std::cout << hessian << std::endl;
    #endif

    // Return Hessian
    return hessian;
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
