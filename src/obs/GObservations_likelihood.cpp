/***************************************************************************
 *         GObservations_likelihood.cpp - Likelihood function class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2016 by Juergen Knoedlseder                         *
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
#include "GFilename.hpp"
#include "GEvent.hpp"
#include "GEventList.hpp"
#include "GEventCube.hpp"
#include "GEventBin.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableStringCol.hpp"
#include "GFitsTableDoubleCol.hpp"

/* __ OpenMP section _____________________________________________________ */
#ifdef _OPENMP
#include <omp.h>
#endif

/* __ Method name definitions ____________________________________________ */
#define G_EVAL             "GObservations::likelihood::eval(GOptimizerPars&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
//#define G_USE_HESSIAN

/* __ Debug definitions __________________________________________________ */
//#define G_EVAL_TIMING     //!< Perform optimizer timing (0=no, 1=yes)
//#define G_EVAL_DEBUG      //!< Perform optimizer debugging (0=no, 1=yes)
//#define G_HESSIAN         //!< Debug Hessian computation

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
        if (m_gradient  != NULL) delete m_gradient;
        if (m_curvature != NULL) delete m_curvature;

        // Initialise value, gradient vector and curvature matrix
        m_value     = 0.0;
        m_npred     = 0.0;
        m_gradient  = new GVector(npars);
        m_curvature = new GMatrixSparse(npars,npars);

        // Set stack size and number of entries
        int stack_size  = (2*npars > 100000) ? 2*npars : 100000;
        int max_entries =  2*npars;
        m_curvature->stack_init(stack_size, max_entries);

        // Allocate vectors to save working variables of each thread
        std::vector<GVector*>       vect_cpy_grad;
        std::vector<GMatrixSparse*> vect_cpy_curvature;
        std::vector<double*>        vect_cpy_value;
        std::vector<double*>        vect_cpy_npred;

        // Here OpenMP will paralellize the execution. The following code will
        // be executed by the differents threads. In order to avoid protecting
        // attributes ( m_value,m_npred, m_gradient and m_curvature), each thread
        // works with its own working variables (cpy_*). When a thread starts,
        // we add working variables in a vector (vect_cpy_*). When computation
        // is finished we just add all elements contain in the vector to the
        // attributes value.
        #pragma omp parallel
        {
            // Allocate and initialize variable copies for multi-threading
            GModels        cpy_model(m_this->models());
            GVector*       cpy_gradient  = new GVector(npars);
            GMatrixSparse* cpy_curvature = new GMatrixSparse(npars,npars);
            double*        cpy_npred     = new double(0.0);
            double*        cpy_value     = new double(0.0);

            // Set stack size and number of entries
            cpy_curvature->stack_init(stack_size, max_entries);

            // Push variable copies into vector. This is a critical zone to
            // avoid multiple thread pushing simultaneously.
            #pragma omp critical
            {
                vect_cpy_grad.push_back(cpy_gradient);
                vect_cpy_curvature.push_back(cpy_curvature); 
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
                                                           cpy_curvature,
                                                           cpy_npred);

            } // endfor: looped over observations

            // Release stack
            cpy_curvature->stack_destroy();

        } // end pragma omp parallel

        // Now the computation is finished, update attributes.
        // For each omp section, a thread will be created.
        #pragma omp sections
        {
            #pragma omp section
            {
                for (int i = 0; i < vect_cpy_curvature.size() ; ++i) {
                    *m_curvature += *(vect_cpy_curvature.at(i));
                    delete vect_cpy_curvature.at(i);
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
        m_curvature->stack_destroy();

    } while(0); // endwhile: main loop

    // Copy over the parameter gradients for all parameters that are
    // free (so that we can access the gradients from outside)
    for (int ipar = 0; ipar < pars.size(); ++ipar) {
        if (pars[ipar]->is_free()) {
            GOptimizerPar* par = const_cast<GOptimizerPar*>(pars[ipar]);
            par->factor_gradient((*m_gradient)[ipar]);
        }
    }

    // Optionally use Hessian instead of curvature matrix
    #if defined(G_USE_HESSIAN)
    *m_curvature = hessian(pars);
    #endif

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


/***********************************************************************//**
 * @brief Compute Hessian matrix
 *
 * @param[in] pars Optimizer parameters.
 *
 * @return Hessian matrix.
 ***************************************************************************/
GMatrixSparse GObservations::likelihood::hessian(const GOptimizerPars& pars)
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
    std::cout << "GObservations::likelihood::hessian: ";
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
    std::cout << "GObservations::likelihood::hessian: " << std::endl;
    std::cout << hessian << std::endl;
    #endif

    // Return Hessian
    return hessian;
}


/***********************************************************************//**
 * @brief Compute covariance matrix
 *
 * @return Covariance matrix.
 *
 * The covariance matrix is calculated as the inverse of the curvature
 * matrix. The method retrieves the model scaling factors from the models
 * that are stored with the observations and multiplies the covariance
 * matrix with the scaling factors of the model parameters. Hence the
 * covariance matrix is covariance with respect to the true model values.
 ***************************************************************************/
GMatrixSparse GObservations::likelihood::covariance(void) const
{
    // Compute covariance matrix (circumvent const correctness)
    GMatrixSparse covmat =
        const_cast<GObservations::likelihood*>(this)->curvature()->invert();

    // Get model parameters
    GOptimizerPars pars = m_this->m_models.pars();

    // Multiply covariance matrix elements with scale factors
    for (int row = 0; row < covmat.rows(); ++row) {
        for (int col = 0; col < covmat.columns(); ++col) {
            covmat(row,col) *= pars[row]->scale() * pars[col]->scale();
        }
    }

    // Return covariance matrix
    return covmat;
}


/***********************************************************************//**
 * @brief Save likelihood fit results into FITS file.
 *
 * @param[in] filename FITS filename.
 *
 * Saves the likelihood fit results into a FITS file. For the moment the
 * method only writes the covariance matrix.
 ***************************************************************************/
void GObservations::likelihood::save(const GFilename& filename) const
{
    // Get covariance matrix
    GMatrixSparse covmat = covariance();

    // Create covariance matrix entry names
    std::vector<std::string> covmat_entries;
    for (int i = 0 ; i < m_this->m_models.size(); ++i) {
	    for (int j = 0; j < m_this->m_models[i]->size(); ++j) {
            covmat_entries.push_back(m_this->m_models[i]->at(j).name() + "(" +
                                     m_this->m_models[i]->name() + ")");
        }
    }

    // Create binary table and columns
    int size = covmat_entries.size();
    GFitsBinTable       covmat_table;
    GFitsTableStringCol par("Parameters", 1, 50, size);
    GFitsTableDoubleCol cov("Covariance", 1, size*size);

    // Fill tables
    int counter = 0;
    for (int i = 0; i < size; ++i) {
        par(0, i) = covmat_entries[i];
        for (int j = 0; j < size; ++j) {
            cov(0, counter) = covmat(i,j);
            ++counter;
        }
    }

    // Set dimension for covariance matrix column
    std::vector<int> dim;
    dim.push_back(size);
    dim.push_back(size);
    cov.dim(dim);

    // Append columns to table
    covmat_table.append(par);
    covmat_table.append(cov);

    // Set extension name
    covmat_table.extname("Covariance Matrix");

    // Allocate FITS object
    GFits fits;

    // Append covariance matrix table to FITS object
    fits.append(covmat_table);

    // Save FITS file to disk
    fits.saveto(filename, true);

    // Close FITS object
    fits.close();

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
    m_curvature = NULL;

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
    m_value = fct.m_value;
    m_npred = fct.m_npred;
    m_this  = fct.m_this;

    // Clone gradient if it exists
    if (fct.m_gradient != NULL) {
        m_gradient = new GVector(*fct.m_gradient);
    }

    // Clone curvature matrix if it exists
    if (fct.m_curvature != NULL) {
        m_curvature = new GMatrixSparse(*fct.m_curvature);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GObservations::likelihood::free_members(void)
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
