/***************************************************************************
 *            GObservation.cpp - Abstract observation base class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2014 by Juergen Knoedlseder                         *
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
 * @file GObservation.cpp
 * @brief Abstract observation base class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GObservation.hpp"
#include "GModelSky.hpp"
#include "GModelData.hpp"
#include "GIntegral.hpp"
#include "GDerivative.hpp"
#include "GTools.hpp"
#include "GEventCube.hpp"
#include "GEventList.hpp"
#include "GEventBin.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_LIKELIHOOD           "GObservation::likelihood(GModels&, GVector*,"\
                                                  " GMatrixSparse*, double*)"
#define G_MODEL                   "GObservation::model(GModels&, GPointing&,"\
                                    " GInstDir&, GEnergy&, GTime&, GVector*)"
#define G_EVENTS                                     "GObservation::events()"
#define G_NPRED_TEMP                 "GObservation::npred_temp(GModel&, int)"
#define G_NPRED_SPEC              "GObservation::npred_spec(GModel&, GTime&)"
#define G_NPRED_SPAT       "GObservation::npred_spat(GModel&, int, GEnergy&,"\
                                                                   " GTime&)"
#define G_NPRED_KERN            "GObservation::npred_kern(GModel&, GSkyDir&,"\
                                             " GEnergy&, GTime&, GPointing&)"
#define G_NPRED_GRAD_TEMP       "GObservation::npred_grad_temp(GModel&, int)"
#define G_NPRED_GRAD_SPEC       "GObservation::npred_grad_spec(GModel&, int,"\
                                                                   " GTime&)"
#define G_NPRED_GRAD_SPAT       "GObservation::npred_grad_spat(GModel&, int,"\
                                                         " GEnergy&, GTime&)"
#define G_NPRED_GRAD_KERN       "GObservation::npred_grad_kern(GModel&, int,"\
                                   " GSkyDir&, GEnergy&, GTime&, GPointing&)"

/* __ Constants __________________________________________________________ */
const double minmod = 1.0e-100;                      //!< Minimum model value
const double minerr = 1.0e-100;                //!< Minimum statistical error

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_LN_ENERGY_INT   //!< ln(E) variable substitution for integration
//#define G_GRAD_RIDDLER  //!< Use Riddler's method for computing derivatives

/* __ Debug definitions __________________________________________________ */
//#define G_OPT_DEBUG


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GObservation::GObservation(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] obs Observation.
 *
 * Instantiate the class by copying from an existing observation.
 ***************************************************************************/
GObservation::GObservation(const GObservation& obs)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(obs);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GObservation::~GObservation(void)
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
 * @param[in] obs Observation.
 * @return Observation.
 *
 * Assign observation.
 ***************************************************************************/
GObservation& GObservation::operator=(const GObservation& obs)
{
    // Execute only if object is not identical
    if (this != &obs) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(obs);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Compute likelihood function
 *
 * @param[in] models Models.
 * @param[in,out] gradient Pointer to gradients.
 * @param[in,out] curvature Pointer to curvature matrix.
 * @param[in,out] npred Pointer to Npred value.
 * @return Likelihood.
 *
 * Computes the likelihood for a specified set of models. The method also
 * returns the gradients, the curvature matrix, and the number of events
 * that are predicted by all models.
 ***************************************************************************/
double GObservation::likelihood(const GModels& models,
                                GVector*       gradient,
                                GMatrixSparse* curvature,
                                double*        npred) const
{
    // Initialise likelihood value
    double value = 0.0;

    // Extract statistics for this observation
    std::string statistics = gammalib::toupper(this->statistics());

    // Unbinned analysis
    if (dynamic_cast<const GEventList*>(events()) != NULL) {

        // Poisson statistics
        if (statistics == "POISSON") {

            // Update the log-likelihood
            value = likelihood_poisson_unbinned(models,
                                                gradient,
                                                curvature,
                                                npred);

        } // endif: Poisson statistics

        // ... otherwise throw an exception
        else {
            throw GException::invalid_statistics(G_LIKELIHOOD, statistics,
                  "Unbinned optimization requires Poisson statistics.");
        }

    } // endif: unbinned analysis

    // ... or binned analysis
    else {

        // Poisson statistics
        if (statistics == "POISSON") {
            value = likelihood_poisson_binned(models,
                                              gradient,
                                              curvature,
                                              npred);
        }

        // ... or Gaussian statistics
        else if (statistics == "GAUSSIAN") {
            value = likelihood_gaussian_binned(models,
                                              gradient,
                                              curvature,
                                              npred);
        }

        // ... or unsupported
        else {
            throw GException::invalid_statistics(G_LIKELIHOOD, statistics,
                  "Binned optimization requires Poisson or Gaussian"
                  " statistics.");
        }

    } // endelse: binned analysis

    // Return likelihood
    return value;
}


/***********************************************************************//**
 * @brief Return model value and (optionally) gradient
 *
 * @param[in] models Model descriptor.
 * @param[in] event Observed event.
 * @param[out] gradient Pointer to gradient vector (optional).
 *
 * @exception GException::gradient_par_mismatch
 *            Dimension of gradient vector mismatches number of parameters.
 *
 * Implements generic model and gradient evaluation for an observation. The
 * model gives the probability for an event to occur with a given instrument
 * direction, at a given energy and at a given time. The gradient is the
 * parameter derivative of this probability. If NULL is passed for the
 * gradient vector, then gradients will not be computed.
 *
 * The method will only operate on models for which the list of instruments
 * and observation identifiers matches those of the observation. Models that
 * do not match will be skipped.
 ***************************************************************************/
double GObservation::model(const GModels& models, const GEvent& event,
                           GVector* gradient) const
{
    // Verify that gradient vector and models have the same dimension
    #if defined(G_RANGE_CHECK)
    if (gradient != NULL) {
        if (models.npars() != gradient->size()) {
            throw GException::gradient_par_mismatch(G_MODEL, 
                                                    gradient->size(),
                                                    models.npars());
        }
    }
    #endif

    // Initialise
    double model = 0.0;    // Reset model value
    int    igrad = 0;      // Reset gradient counter

    // If gradient is available then reset gradient vector elements to 0
    if (gradient != NULL) {
        (*gradient) = 0.0;
    }

    // Loop over models
    for (int i = 0; i < models.size(); ++i) {

        // Get model pointer. Continue only if pointer is valid
        const GModel* mptr = models[i];
        if (mptr != NULL) {

            // Continue only if model applies to specific instrument and
            // observation identifier
            if (mptr->is_valid(instrument(), id())) {

                // Compute value and add to model. If energy dispersion
                // is used, don't compute model gradients as we cannot
                // use them. This is somehow a kluge, but makes the
                // code faster
                if (response().use_edisp()) {
                    model += mptr->eval(event, *this);
                }
                else {
                    model += mptr->eval_gradients(event, *this);
                }

                // Optionally determine model gradients
                if (gradient != NULL) {
                    for (int k = 0; k < mptr->size(); ++k) {
                        (*gradient)[igrad+k] = model_grad(*mptr, event, k);
                    }
                }

            } // endif: model component was valid for instrument

            // Increment parameter counter for gradients
            igrad += mptr->size();

        } // endif: model was valid

    } // endfor: Looped over models

    // Return
    return model;
}


/***********************************************************************//**
 * @brief Return total number (and optionally gradient) of predicted counts
 *        for all models
 *
 * @param[in] models Models.
 * @param[out] gradient Model parameter gradients (optional).
 *
 * @exception GException::gradient_par_mismatch
 *            Dimension of gradient vector mismatches number of parameters.
 *
 * Returns the total number of predicted counts within the analysis region.
 * If NULL is passed for the gradient vector then gradients will not be
 * computed.
 *
 * The method will only operate on models for which the list of instruments
 * and observation identifiers matches those of the observation. Models that
 * do not match will be skipped.
 ***************************************************************************/
double GObservation::npred(const GModels& models, GVector* gradient) const
{
    // Verify that gradient vector and models have the same dimension
    #if defined(G_RANGE_CHECK)
    if (gradient != NULL) {
        if (models.npars() != gradient->size()) {
            throw GException::gradient_par_mismatch(G_MODEL, 
                                                    gradient->size(),
                                                    models.npars());
        }
    }
    #endif

    // Initialise
    double npred = 0.0;    // Reset predicted number of counts
    int    igrad = 0;      // Reset gradient counter

    // If gradient is available then reset gradient vector elements to 0
    if (gradient != NULL) {
        (*gradient) = 0.0;
    }

    // Loop over models
    for (int i = 0; i < models.size(); ++i) {

        // Get model pointer. Continue only if pointer is valid
        const GModel* mptr = models[i];
        if (mptr != NULL) {

            // Continue only if model applies to specific instrument and
            // observation identifier
            if (mptr->is_valid(instrument(), id())) {

                // Determine Npred for model
                npred += npred_temp(*mptr);

                // Optionally determine Npred gradients
                if (gradient != NULL) {
                    for (int k = 0; k < mptr->size(); ++k) {
                        (*gradient)[igrad+k] = npred_grad(*mptr, k);
                    }
                }

            } // endif: model component was valid for instrument

            // Increment parameter counter for gradient
            igrad += mptr->size();

        } // endif: model was valid

    } // endfor: Looped over models

    // Return prediction
    return npred;
}


/***********************************************************************//**
 * @brief Returns parameter gradient of model for a given event
 *
 * @param[in] model Model.
 * @param[in] event Event.
 * @param[in] ipar Parameter index for which gradient should be returned.
 *
 * This method uses a robust but dumb method to estimate parameter
 * gradients that have not been provided by the model. We use here a dumb
 * method as this method is likely used for spatial model parameters, and
 * the spatial model may eventually be noisy due to numerical integration
 * limits.
 *
 * The step size for the dumb method has been fixed to 0.0002, which
 * corresponds to abound 1 arcsec for parameters that are given in degrees.
 * The reasoning behind this value is that parameters that use numerical
 * gradients are typically angles, such as for example the position, and
 * we want to achieve arcsec precision with this method.
 *
 * @todo Implement a more precise numerical derivation scheme. This needs
 *       a deep investigation as the scheme needs to be precise but also
 *       robust (not sensitive to noise in the function).
 * @todo For the Riddler method, we simply remove any parameter boundaries
 *       here for the computation to avoid any out of boundary errors.
 *       We may have models, however, for which out of bound parameters lead
 *       to illegal computations, such as division by zero or taking the
 *       square root of negative values.
 *       I cannot see any elegant method to catch this at this level.
 *       Eventually, the higher level method should avoid going in a
 *       parameter domain that is not defined. 
 ***************************************************************************/
double GObservation::model_grad(const GModel& model,
                                const GEvent& event,
                                const int&    ipar) const
{
    // Initialise gradient
    double grad = 0.0;

    // Get reference on model parameter to avoid repreated operator access
    const GModelPar& par = model[ipar];

    // Compute gradient only if parameter is free
    if (par.is_free()) {

        // If model has a gradient then use it, unless we have energy
        // dispersion. For energy dispersion, no gradients are available
        // as we have not implemented code that integrates the gradients
        // over the energy dispersion. Here it's simpler to just use
        // numerical gradients.
        if (par.has_grad() && !response().use_edisp()) {
            grad = par.factor_gradient();
        }

        // ... otherwise compute it numerically
        else {

            // Get non-const model pointer
            GModelPar* ptr = const_cast<GModelPar*>(&par);

            // Save current model parameter
            GModelPar current = par;

            // Get actual parameter value
            double x = par.factor_value();

            // Set fixed step size for computation of derivative.
            // By default, the step size is fixed to 0.0002, but if this would
            // violate a boundary, dx is reduced accordingly. In case that x
            // is right on the boundary, x is displaced slightly from the
            // boundary to allow evaluation of the derivative.
            #if !defined(G_GRAD_RIDDLER)
            const double step_size = 0.0002;
            double       dx        = step_size;
            if (model[ipar].has_min()) {
                double dx_min = x - par.factor_min();
                if (dx_min == 0.0) {
                    dx = step_size * x;
                    if (dx == 0.0) {
                        dx = step_size;
                    }
                    x += dx; 
                }
                else if (dx_min < dx) {
                    dx = dx_min;
                }
            }
            if (model[ipar].has_max()) {
                double dx_max = par.factor_max() - x;
                if (dx_max == 0.0) {
                    dx = step_size * x;
                    if (dx == 0.0) {
                        dx = step_size;
                    }
                    x -= dx; 
                }
                else if (dx_max < dx) {
                    dx = dx_max;
                }
            }
            #endif

            // Remove any boundaries to avoid limitations
            ptr->remove_range();

            // Setup derivative function
            GObservation::model_func function(this, model, event, ipar);

            // Get derivative. We use a fixed step size here that has been
            // checked on spatial parameters of models
            GDerivative derivative(&function);
            #if defined(G_GRAD_RIDDLER)
            grad = derivative.value(x);
            #else
            grad = derivative.difference(x, dx);
            #endif

            // Restore current model parameter
            *ptr = current;

        } // endelse: computed gradient numerically

    } // endif: model parameter was free

    // Return gradient
    return grad;
}


/***********************************************************************//**
 * @brief Returns parameter gradient of Npred
 *
 * @param[in] model Gamma-ray source model.
 * @param[in] ipar Parameter index for which gradient should be returned.
 *
 * Computes
 * \f[\frac{{\rm d} N_{\rm pred}}{{\rm d} a_i}\f]
 * where
 * \f[N_{\rm pred} = \int_{\rm GTI} \int_{E_{\rm bounds}} \int_{\rm ROI}
 *    S(\vec{p}, E, t) PSF(\vec{p'}, E', t' | \vec{d}, \vec{p}, E, t) \,
 *    {\rm d}\vec{p'} {\rm d}E' {\rm d}t'\f]
 * and
 * \f$a_i\f$ is the model parameter \f$i\f$.
 * Furthermore,
 * \f$S(\vec{p}, E, t)\f$ is the source model,
 * \f$PSF(\vec{p'}, E', t' | \vec{d}, \vec{p}, E, t)\f$ is the point
 * spread function,
 * \f$\vec{p'}\f$ is the measured photon direction,
 * \f$E'\f$ is the measured photon energy,
 * \f$t'\f$ is the measured photon arrival time,
 * \f$\vec{p}\f$ is the true photon arrival direction,
 * \f$E\f$ is the true photon energy,
 * \f$t\f$ is the true photon arrival time, and
 * \f$d\f$ is the instrument pointing.
 *
 * This method uses a robust but dumb method to estimate gradients. This
 * method has turned out more robust then the Riddler's method implement
 * by the GDerivative::value() method.
 *
 * The step size for the dumb method has been fixed to 0.0002, which
 * corresponds to abound 1 arcsec for parameters that are given in degrees.
 * The reasoning behind this value is that parameters that use numerical
 * gradients are typically angles, such as for example the position, and
 * we want to achieve arcsec precision with this method.
 *
 * @todo Implement a more precise numerical derivation scheme. This needs
 *       a deep investigation as the scheme needs to be precise but also
 *       robust (not sensitive to noise in the function).
 * @todo For Riddler's method we simply remove any parameter boundaries here
 *       for the computation to avoid any out of boundary errors. We may have
 *       models, however, for which out of bound parameters lead to illegal
 *       computations, such as division by zero or taking the square root of
 *       negative values.
 *       I cannot see any elegant method to catch this at this level.
 *       Eventually, the higher level method should avoid going in a
 *       parameter domain that is not defined. 
 ***************************************************************************/
double GObservation::npred_grad(const GModel& model, const int& ipar) const
{
    // Initialise result
    double grad = 0.0;

    // Compute gradient only if parameter is free
    if (model[ipar].is_free()) {

        // Get non-const model pointer (circumvent const correctness)
        GModel* ptr = const_cast<GModel*>(&model);

        // Save current model parameter
        GModelPar current = (*ptr)[ipar];

        // Get actual parameter value
        double x = model[ipar].factor_value();

        // Determine fixed step size for computation of derivative.
        // By default, the step size is fixed to 0.0002, but if this would
        // violate a boundary, dx is reduced accordingly. In case that x
        // is right on the boundary, x is displaced slightly from the
        // boundary to allow evaluation of the derivative.
        #if !defined(G_GRAD_RIDDLER)
        const double step_size = 0.0002;
        double       dx        = step_size;
        if (model[ipar].has_min()) {
            double dx_min = x - model[ipar].factor_min();
            if (dx_min == 0.0) {
                dx = step_size * x;
                if (dx == 0.0) {
                    dx = step_size;
                }
                x += dx; 
            }
            else if (dx_min < dx) {
                dx = dx_min;
            }
        }
        if (model[ipar].has_max()) {
            double dx_max = model[ipar].factor_max() - x;
            if (dx_max == 0.0) {
                dx = step_size * x;
                if (dx == 0.0) {
                    dx = step_size;
                }
                x -= dx; 
            }
            else if (dx_max < dx) {
                dx = dx_max;
            }
        }
        #endif

        // Remove any boundaries to avoid limitations
        (*ptr)[ipar].remove_range();

        // Setup derivative function
        GObservation::npred_func function(this, model, ipar);

        // Get derivative.
        GDerivative derivative(&function);
        #if defined(G_GRAD_RIDDLER)
        grad = derivative.value(x);
        #else
        grad = derivative.difference(x, dx);
        #endif

        // Restore current model parameter
        (*ptr)[ipar] = current;

    } // endif: model parameter was free

    // Return result
    return grad;
}


/***********************************************************************//**
 * @brief Set event container
 *
 * @param[in] events Event container.
 *
 * Set the event container for this observation by cloning the container
 * specified in the argument.
 ***************************************************************************/
void GObservation::events(const GEvents& events)
{
    // Remove an existing event container
    if (m_events != NULL) delete m_events;

    // Clone events
    m_events = events.clone();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return event container
 *
 * @exception GException::no_events
 *            No event container defined for observation.
 *
 * Returns pointer to event container.
 ***************************************************************************/
const GEvents* GObservation::events(void) const
{
    // Throw an exception if the event container is not valid
    if (m_events == NULL) {
        throw GException::no_events(G_EVENTS);
    }

    // Return pointer to event container
    return m_events;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GObservation::init_members(void)
{
    // Initialise members
    m_name.clear();
    m_id.clear();
    m_statistics = "Poisson";
    m_events     = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] obs Observation.
 *
 * Copy members from an observation.
 ***************************************************************************/
void GObservation::copy_members(const GObservation& obs)
{
    // Copy members
    m_name       = obs.m_name;
    m_id         = obs.m_id;
    m_statistics = obs.m_statistics;

    // Clone members
    m_events = (obs.m_events != NULL) ? obs.m_events->clone() : NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GObservation::free_members(void)
{
    // Free members
    if (m_events != NULL) delete m_events;

    // Signal free pointers
    m_events = NULL;

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                           Likelihood methods                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Evaluate log-likelihood function for Poisson statistics and
 *        unbinned analysis (version with working arrays)
 *
 * @param[in] models Models.
 * @param[in,out] gradient Gradient.
 * @param[in,out] curvature Curvature matrix.
 * @param[in,out] npred Number of predicted events.
 * @return Likelihood value.
 *
 * This method evaluates the -(log-likelihood) function for parameter
 * optimisation using unbinned analysis and Poisson statistics.
 * The -(log-likelihood) function is given by
 *
 * \f$L = N_{\rm pred} - \sum_i \log e_i\f$
 *
 * where
 * \f$N_{\rm pred}\f$ is the number of events predicted by the model, and
 * the sum is taken over all events. This method also computes the
 * parameter gradients
 * \f$\delta L/dp\f$
 * and the curvature matrix
 * \f$\delta^2 L/dp_1 dp_2\f$.
 ***************************************************************************/
double GObservation::likelihood_poisson_unbinned(const GModels& models,
                                                 GVector*       gradient,
                                                 GMatrixSparse* curvature,
                                                 double*        npred) const
{
    // Initialise likelihood value
    double value = 0.0;

    // Get number of parameters
    int npars = gradient->size();

    // Allocate some working arrays
    int*    inx    = new int[npars];
    double* values = new double[npars];
    GVector wrk_grad(npars);

    // Determine Npred value and gradient for this observation
    double npred_value = this->npred(models, &wrk_grad);

    // Update likelihood, Npred and gradient
    value     += npred_value;
    *npred    += npred_value;
    *gradient += wrk_grad;

    // Iterate over all events
    for (int i = 0; i < events()->size(); ++i) {

        // Get event pointer
        const GEvent* event = (*events())[i];

        // Get model and derivative
        double model = this->model(models, *event, &wrk_grad);

        // Skip bin if model is too small (avoids -Inf or NaN gradients)
        if (model <= minmod) {
            continue;
        }

        // Create index array of non-zero derivatives and initialise working
        // array
        int ndev = 0;
        for (int i = 0; i < npars; ++i) {
            values[i] = 0.0;
            if (wrk_grad[i] != 0.0 && !gammalib::is_infinite(wrk_grad[i])) {
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
            double       g       = wrk_grad[jpar];
            double       fa_i    = fa * g;

            // Update gradient.
            (*gradient)[jpar] -= fb * g;

            // Loop over rows
            register int* ipar = inx;

            for (register int idev = 0; idev < ndev; ++idev, ++ipar) {
                values[idev] = fa_i * wrk_grad[*ipar];
            }

            // Add column to matrix
            curvature->add_to_column(jpar, values, inx, ndev);

        } // endfor: looped over columns

    } // endfor: iterated over all events

    // Free temporary memory
    if (values != NULL) delete [] values;
    if (inx    != NULL) delete [] inx;

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate log-likelihood function for Poisson statistics and
 *        binned analysis (version with working arrays)
 *
 * @param[in] models Models.
 * @param[in,out] gradient Gradient.
 * @param[in,out] curvature Curvature matrix.
 * @param[in,out] npred Number of predicted events.
 * @return Likelihood value.
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
double GObservation::likelihood_poisson_binned(const GModels& models,
                                             GVector*       gradient,
                                             GMatrixSparse* curvature,
                                             double*        npred) const
{
    // Initialise likelihood value
    double value = 0.0;

    // Initialise statistics
    #if defined(G_OPT_DEBUG)
    int    n_bins        = 0;
    int    n_used        = 0;
    int    n_small_model = 0;
    int    n_zero_data   = 0;
    double sum_data      = 0.0;
    double sum_model     = 0.0;
    double init_value    = value;
    #endif

    // Get number of parameters
    int npars = gradient->size();

    // Allocate some working arrays
    int*    inx    = new int[npars];
    double* values = new double[npars];
    GVector wrk_grad(npars);

    // Iterate over all bins
    for (int i = 0; i < events()->size(); ++i) {

        // Update number of bins
        #if defined(G_OPT_DEBUG)
        n_bins++;
        #endif

        // Get event pointer
        const GEventBin* bin =
            (*(static_cast<GEventCube*>(const_cast<GEvents*>(events()))))[i];

        // Get number of counts in bin
        double data = bin->counts();

        // Get model and derivative
        double model = this->model(models, *bin, &wrk_grad);

        // Multiply model by bin size
        model *= bin->size();

        // Skip bin if model is too small (avoids -Inf or NaN gradients)
        if (model <= minmod) {
            #if defined(G_OPT_DEBUG)
            n_small_model++;
            #endif
            continue;
        }

        // Update statistics
        #if defined(G_OPT_DEBUG)
        n_used++;
        sum_data  += data;
        sum_model += model;
        #endif

        // Update Npred
        *npred += model;

        // Multiply gradient by bin size
        wrk_grad *= bin->size();

        // Create index array of non-zero derivatives and initialise working
        // array
        int ndev = 0;
        for (int i = 0; i < npars; ++i) {
            values[i] = 0.0;
            if (wrk_grad[i] != 0.0 && !gammalib::is_infinite(wrk_grad[i])) {
                inx[ndev] = i;
                ndev++;
            }
        }

        // Update gradient vector and curvature matrix. To avoid
        // unneccessary computations we distinguish the case where
        // data>0 and data=0. The second case requires much less
        // computation since it does not contribute to the curvature
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
                double       g       = wrk_grad[jpar];
                double       fa_i    = fa * g;

                // Update gradient
                (*gradient)[jpar] += fc * g;

                // Loop over rows
                register int* ipar = inx;
                for (register int idev = 0; idev < ndev; ++idev, ++ipar) {
                    values[idev] = fa_i * wrk_grad[*ipar];
                }

                // Add column to matrix
                curvature->add_to_column(jpar, values, inx, ndev);

            } // endfor: looped over columns

        } // endif: data was > 0

        // ... handle now data=0
        else {

            // Update statistics
            #if defined(G_OPT_DEBUG)
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
                (*gradient)[*ipar] += wrk_grad[*ipar];
            }

        } // endif: data was 0

    } // endfor: iterated over all events

    // Free temporary memory
    if (values != NULL) delete [] values;
    if (inx    != NULL) delete [] inx;

    // Dump statistics
    #if defined(G_OPT_DEBUG)
    std::cout << "Number of bins: " << n_bins << std::endl;
    std::cout << "Number of bins used for computation: " << n_used << std::endl;
    std::cout << "Number of bins excluded due to small model: " << n_small_model << std::endl;
    std::cout << "Number of bins with zero data: " << n_zero_data << std::endl;
    std::cout << "Sum of data: " << sum_data << std::endl;
    std::cout << "Sum of model: " << sum_model << std::endl;
    std::cout << "Initial statistics: " << init_value << std::endl;
    std::cout << "Statistics: " << m_value-init_value << std::endl;
    #endif

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate log-likelihood function for Gaussian statistics and
 *        binned analysis (version with working arrays)
 *
 * @param[in] models Models.
 * @param[in,out] gradient Gradient.
 * @param[in,out] curvature Curvature matrix.
 * @param[in,out] npred Number of predicted events.
 * @return Likelihood value.
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
double GObservation::likelihood_gaussian_binned(const GModels& models,
                                                GVector*       gradient,
                                                GMatrixSparse* curvature,
                                                double*        npred) const
{
    // Initialise likelihood value
    double value = 0.0;

    // Get number of parameters
    int npars = gradient->size();

    // Allocate some working arrays
    int*    inx    = new int[npars];
    double* values = new double[npars];
    GVector wrk_grad(npars);

    // Iterate over all bins
    for (int i = 0; i < events()->size(); ++i) {

        // Get event pointer
        const GEventBin* bin =
            (*(static_cast<GEventCube*>(const_cast<GEvents*>(events()))))[i];

        // Get number of counts in bin
        double data = bin->counts();

        // Get statistical uncertainty
        double sigma = bin->error();

        // Skip bin if statistical uncertainty is too small
        if (sigma <= minerr) {
            continue;
        }

        // Get model and derivative
        double model = this->model(models, *bin, &wrk_grad);

        // Multiply model by bin size
        model *= bin->size();

        // Skip bin if model is too small (avoids -Inf or NaN gradients)
        if (model <= minmod) {
            continue;
        }

        // Update Npred
        *npred += model;

        // Multiply gradient by bin size
        wrk_grad *= bin->size();

        // Create index array of non-zero derivatives and initialise working
        // array
        int ndev = 0;
        for (int i = 0; i < npars; ++i) {
            values[i] = 0.0;
            if (wrk_grad[i] != 0.0 && !gammalib::is_infinite(wrk_grad[i])) {
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
            double       fa_i = wrk_grad[jpar] * weight;

            // Update gradient
            (*gradient)[jpar] -= fa * fa_i;

            // Loop over rows
            register int* ipar = inx;
            for (register int idev = 0; idev < ndev; ++idev, ++ipar) {
                values[idev] = fa_i * wrk_grad[*ipar];
            }

            // Add column to matrix
            curvature->add_to_column(jpar, values, inx, ndev);

        } // endfor: looped over columns

    } // endfor: iterated over all events

    // Free temporary memory
    if (values != NULL) delete [] values;
    if (inx    != NULL) delete [] inx;

    // Return
    return value;
}


/*==========================================================================
 =                                                                         =
 =                         Model gradient methods                          =
 =                                                                         =
 ==========================================================================*/



/***********************************************************************//**
 * @brief Model function evaluation for gradient computation
 *
 * @param[in] x Function value.
 ***************************************************************************/
double GObservation::model_func::eval(const double& x)
{
    // Get non-const model pointer (circumvent const correctness)
    GModel* model = const_cast<GModel*>(m_model);

    // Set value
    (*model)[m_ipar].factor_value(x);

    // Compute model value
    double value = model->eval(*m_event, *m_parent);

    // Return value
    return value;
}


/*==========================================================================
 =                                                                         =
 =                         Npred computation methods                       =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Npred function evaluation for gradient computation
 *
 * @param[in] x Function value.
 ***************************************************************************/
double GObservation::npred_func::eval(const double& x)
{
    // Get non-const model pointer (circumvent const correctness)
    GModel* model = const_cast<GModel*>(m_model);

    // Set value
    (*model)[m_ipar].factor_value(x);

    // Compute Npred value
    double npred = m_parent->npred_temp(*model);

    // Return value
    return npred;
}


/***********************************************************************//**
 * @brief Temporally integrates spatially and spectrally integrated Npred
 *        kernel
 *
 * @param[in] model Gamma-ray source model.
 *
 * @exception GException::gti_invalid
 *            Good Time Interval is invalid.
 *
 * Computes
 * \f[N_{\rm pred} = \int_{\rm GTI} \int_{E_{\rm bounds}} \int_{\rm ROI}
 *    S(\vec{p}, E, t) PSF(\vec{p'}, E', t' | \vec{d}, \vec{p}, E, t) \,
 *    {\rm d}\vec{p'} {\rm d}E' {\rm d}t'\f]
 * where
 * \f$S(\vec{p}, E, t)\f$ is the source model,
 * \f$PSF(\vec{p'}, E', t' | \vec{d}, \vec{p}, E, t)\f$ is the point
 * spread function,
 * \f$\vec{p'}\f$ is the measured photon direction,
 * \f$E'\f$ is the measured photon energy,
 * \f$t'\f$ is the measured photon arrival time,
 * \f$\vec{p}\f$ is the true photon arrival direction,
 * \f$E\f$ is the true photon energy,
 * \f$t\f$ is the true photon arrival time, and
 * \f$d\f$ is the instrument pointing.
 *
 * \f${\rm GTI}\f$ are the Good Time Intervals that are stored in the
 * GObservation::m_gti member.
 * Note that MET is used for the time integration interval. This, however,
 * is no specialisation since npred_grad_kern_spec::eval() converts the
 * argument back in a GTime object by assuming that the argument is in MET,
 * hence the correct time system will be used at the end by the method.
 ***************************************************************************/
double GObservation::npred_temp(const GModel& model) const
{
    // Initialise result
    double result = 0.0;


    // Case A: If the model is constant then integrate analytically
    if (model.is_constant()) {

        // Evaluate model at first start time and multiply by ontime
        double ontime = events()->gti().ontime();

        // Integrate only if ontime is positive
        if (ontime > 0.0) {

            // Integration is a simple multiplication by the time
            result = npred_spec(model, events()->gti().tstart()) * ontime;

        }

    } // endif: model was constant

    // ... otherwise integrate temporally
    else {

        // Loop over GTIs
        for (int i = 0; i < events()->gti().size(); ++i) {

            // Set integration interval in seconds
            double tstart = events()->gti().tstart(i).secs();
            double tstop  = events()->gti().tstop(i).secs();

            // Throw exception if time interval is not valid
            if (tstop <= tstart) {
                throw GException::gti_invalid(G_NPRED_TEMP,
                                              events()->gti().tstart(i),
                                              events()->gti().tstop(i));
            }

            // Setup integration function
            GObservation::npred_temp_kern integrand(this, &model);
            GIntegral                     integral(&integrand);

            // Do Romberg integration
            result += integral.romb(tstart, tstop);

        } // endfor: looped over GTIs

    } // endelse: integrated temporally

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Integration kernel for npred_temp() method
 *
 * @param[in] x Function value.
 ***************************************************************************/
double GObservation::npred_temp_kern::eval(const double& x)
{
    // Convert argument in native reference in seconds
    GTime time;
    time.secs(x);

    // Return value
    return (m_parent->npred_spec(*m_model, time));
}


/***********************************************************************//**
 * @brief Integrates spatially integrated Npred kernel spectrally
 *
 * @param[in] model Gamma-ray source model.
 * @param[in] obsTime Measured photon arrival time.
 *
 * @exception GException::erange_invalid
 *            Energy range is invalid.
 *
 * Computes
 * \f[N'_{\rm pred} = \int_{E_{\rm bounds}} \int_{\rm ROI}
 *    S(\vec{p}, E, t) PSF(\vec{p'}, E', t' | \vec{d}, \vec{p}, E, t) \,
 *    {\rm d}\vec{p'} {\rm d}E'\f]
 * where
 * \f$S(\vec{p}, E, t)\f$ is the source model,
 * \f$PSF(\vec{p'}, E', t' | \vec{d}, \vec{p}, E, t)\f$ is the point
 * spread function,
 * \f$\vec{p'}\f$ is the measured photon direction,
 * \f$E'\f$ is the measured photon energy,
 * \f$t'\f$ is the measured photon arrival time,
 * \f$\vec{p}\f$ is the true photon arrival direction,
 * \f$E\f$ is the true photon energy,
 * \f$t\f$ is the true photon arrival time, and
 * \f$d\f$ is the instrument pointing.
 *
 * \f$E_{\rm bounds}\f$ are the energy boundaries that are stored in the
 * GObservation::m_ebounds member.
 *
 * @todo Loop also over energy boundaries (is more general; there is no
 *       reason for not doing it).
 ***************************************************************************/
double GObservation::npred_spec(const GModel& model,
                                const GTime&  obsTime) const
{
    // Set integration energy interval in MeV
    GEnergy e_min = events()->ebounds().emin();
    GEnergy e_max = events()->ebounds().emax();

    // If model is sky model and observation uses energy dispersion then
    // add margin
    if (dynamic_cast<const GModelSky*>(&model) != NULL && response().use_edisp()) {
        e_min = response().ebounds_src(e_min).emin();
        e_max = response().ebounds_src(e_max).emax();
    }

    // Get energy interval in MeV
    double emin = e_min.MeV();
    double emax = e_max.MeV();

    // Throw exception if energy range is not valid
    if (emax <= emin) {
        throw GException::erange_invalid(G_NPRED_SPEC, emin, emax);
    }

    // Setup integration function
    GObservation::npred_spec_kern integrand(this, &model, &obsTime);
    GIntegral                     integral(&integrand);

    // Set integration precision
    integral.eps(1.0e-6); // Needed for fluctuating bgd. model !!!

    // Do Romberg integration
    #if defined(G_LN_ENERGY_INT)
    emin = log(emin);
    emax = log(emax);
    #endif
    double result = integral.romb(emin, emax);

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(result) || gammalib::is_infinite(result)) {
        std::cout << "*** ERROR: GObservation::npred_spec:";
        std::cout << " NaN/Inf encountered";
        std::cout << " (result=" << result;
        std::cout << ", emin=" << emin;
        std::cout << ", emax=" << emax;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Integration kernel for npred_spec() method
 *
 * @param[in] x Function value.
 *
 * This method implements the integration kernel needed for the npred_spec()
 * method. If G_LN_ENERGY_INT is defined the energy integration is done
 * logarithmically, i.e. @p x is given in ln(energy) instead of energy.
 ***************************************************************************/
double GObservation::npred_spec_kern::eval(const double& x)
{
    // Set energy
    GEnergy eng;
    #if defined(G_LN_ENERGY_INT)
    double expx = std::exp(x);
    eng.MeV(expx);
    #else
    eng.MeV(x);
    #endif

    // Get function value
    double value = m_model->npred(eng, *m_time, *m_parent);

    // Save value if needed
    #if defined(G_NAN_CHECK)
    double value_out = value;
    #endif

    // Correct for variable substitution
    #if defined(G_LN_ENERGY_INT)
    value *= expx;
    #endif

    // Compile option: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GObservation::npred_spec_kern::eval";
        std::cout << "(x=" << x << "): ";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << " (value_out=" << value_out;
        #if defined(G_LN_ENERGY_INT)
        std::cout << " exp(x)=" << expx;
        #endif
        std::cout << ")" << std::endl;
    }
    #endif

    // Return value
    return value;
}
