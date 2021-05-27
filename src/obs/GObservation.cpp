/***************************************************************************
 *            GObservation.cpp - Abstract observation base class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2021 by Juergen Knoedlseder                         *
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
#include "GTools.hpp"
#include "GException.hpp"
#include "GIntegral.hpp"
#include "GDerivative.hpp"
#include "GVector.hpp"
#include "GMatrixSparse.hpp"
#include "GObservation.hpp"
#include "GResponse.hpp"
#include "GEventCube.hpp"
#include "GEventList.hpp"
#include "GEventBin.hpp"
#include "GModel.hpp"
#include "GModels.hpp"
#include "GModelPar.hpp"
#include "GModelSky.hpp"
#include "GModelData.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_LIKELIHOOD           "GObservation::likelihood(GModels&, GVector*,"\
                                                  " GMatrixSparse*, double*)"
#define G_MODEL1           "GObservation::model(GModels&, GEvent&, GVector*)"
#define G_MODEL2              "GObservation::model(GModels&, GMatrixSparse*)"
#define G_EVENTS                                     "GObservation::events()"
#define G_NPRED                                "GObservation::npred(GModel&)"
#define G_NPRED_SPEC              "GObservation::npred_spec(GModel&, GTime&)"

/* __ Constants __________________________________________________________ */
const double minmod = 1.0e-100;                      //!< Minimum model value
const double minerr = 1.0e-100;                //!< Minimum statistical error

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_LN_ENERGY_INT   //!< ln(E) variable substitution for integration

/* __ Debug definitions __________________________________________________ */
//#define G_OPT_DEBUG                       //!< Debug likelihood computation
//#define G_DEBUG_VECTOR_MODEL              //!< Debug vector model


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
 * @brief Return events
 *
 * @exception GException::no_events
 *            No events allocated for observation.
 *
 * Returns pointer to events.
 ***************************************************************************/
GEvents* GObservation::events(void)
{
    // Throw an exception if the event container is not valid
    if (m_events == NULL) {
        std::string msg = "No events allocated for observation.";
        throw GException::invalid_value(G_EVENTS, msg);
    }

    // Return pointer to events
    return (m_events);
}


/***********************************************************************//**
 * @brief Return events (const version)
 *
 * @exception GException::no_events
 *            No events allocated for observation.
 *
 * Returns const pointer to events.
 ***************************************************************************/
const GEvents* GObservation::events(void) const
{
    // Throw an exception if the event container is not valid
    if (m_events == NULL) {
        std::string msg = "No events allocated for observation.";
        throw GException::invalid_value(G_EVENTS, msg);
    }

    // Return pointer to events
    return (m_events);
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
    // Free existing event container only if it differs from current event
    // container. This prevents unintential deallocation of the argument
    if ((m_events != NULL) && (m_events != &events)) {
        delete m_events;
    }

    // Clone events
    m_events = events.clone();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute likelihood function
 *
 * @param[in] models Models.
 * @param[in,out] gradients Pointer to gradients.
 * @param[in,out] curvature Pointer to curvature matrix.
 * @param[in,out] npred Pointer to Npred value.
 * @return Likelihood.
 *
 * Computes the likelihood for a specified set of models. The method also
 * returns the gradients, the curvature matrix, and the number of events
 * that are predicted by all models.
 ***************************************************************************/
double GObservation::likelihood(const GModels& models,
                                GVector*       gradients,
                                GMatrixSparse* curvature,
                                double*        npred) const
{
    // Initialise likelihood value
    double value = 0.0;

    // Extract statistic for this observation
    std::string statistic = gammalib::toupper(this->statistic());

    // Unbinned analysis
    if (dynamic_cast<const GEventList*>(events()) != NULL) {

        // Poisson statistic
        if ((statistic == "POISSON") || (statistic == "CSTAT")) {

            // Update the log-likelihood
            value = likelihood_poisson_unbinned(models,
                                                gradients,
                                                curvature,
                                                npred);

        } // endif: Poisson statistic

        // ... otherwise throw an exception
        else {
            std::string msg = "Invalid statistic \""+statistic+"\". Unbinned "
                              "optimization requires \"POISSON\" or \"CSTAT\" "
                              "statistic.";
            throw GException::invalid_value(G_LIKELIHOOD, msg);
        }

    } // endif: unbinned analysis

    // ... or binned analysis
    else {

        // Poisson statistic
        if ((statistic == "POISSON") || (statistic == "CSTAT")) {
            value = likelihood_poisson_binned(models,
                                              gradients,
                                              curvature,
                                              npred);
        }

        // ... or Gaussian statistic
        else if ((statistic == "GAUSSIAN")  || (statistic == "CHI2")) {
            value = likelihood_gaussian_binned(models,
                                               gradients,
                                               curvature,
                                               npred);
        }

        // ... or unsupported
        else {
            std::string msg = "Invalid statistic \""+statistic+"\". Binned "
                              "optimization requires \"POISSON\", \"CSTAT\", "
                              "\"GAUSSIAN\" or \"CHI2\" statistic.";
            throw GException::invalid_value(G_LIKELIHOOD, msg);
        }

    } // endelse: binned analysis

    // Return likelihood
    return value;
}


/***********************************************************************//**
 * @brief Return model value and (optionally) gradients
 *
 * @param[in] models Model container.
 * @param[in] event Observed event.
 * @param[out] gradients Pointer to gradient vector (optional).
 * @return Model value.
 *
 * @exception GException::invalid_value
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
double GObservation::model(const GModels& models,
                           const GEvent&  event,
                           GVector*       gradients) const
{
    // Initialise method variables
    double model     = 0.0;    // Reset model value
    int    igrad     = 0;      // Reset gradient counter
    int    grad_size = 0;      // Reset gradient size

    // Set gradient usage flag
    bool use_grad = (gradients != NULL);

    // If gradient is available then reset gradient vector elements to 0
    // and determine vector size
    if (use_grad) {

        // Reset gradient vector
        (*gradients) = 0.0;
        grad_size    = gradients->size();

        // Initialise stack of parameters with gradients
        m_pars_with_gradients.clear();
        m_pars_with_gradients.reserve(models.size());

    } // endif: gradient vector available

    // Loop over models
    for (int i = 0; i < models.size(); ++i) {

        // Get model pointer. Continue only if pointer is valid
        const GModel* mptr = models[i];
        if (mptr != NULL) {

            // Continue only if model applies to specific instrument and
            // observation identifier
            if (mptr->is_valid(instrument(), id())) {

                // Compute value and add to model
                model += mptr->eval(event, *this, use_grad);

                // Optionally determine model gradients
                if (use_grad) {

                    // Make sure that we have a slot for the gradient
                    #if defined(G_RANGE_CHECK)
                    if (igrad+mptr->size() > grad_size) {
                        std::string msg = "Vector has not enough elements to "
                            "store the model parameter gradients. "+
                            gammalib::str(models.npars())+" elements "
                            "requested while vector only contains "+
                            gammalib::str(gradients->size())+" elements.";
                        throw GException::invalid_value(G_MODEL1, msg);
                    }
                    #endif

                    // If model provides the parameter indices that were
                    // updated by the eval() method then use them as this is
                    // less time consuming than looping over all indices
                    if (mptr->has_eval_indices()) {

                        // Get number of relevant parameters
                        int npars = mptr->eval_indices().size();

                        // Loop over all relevant parameters
                        for (int index = 0; index < npars; ++index) {

                            // Retrieve parameter index
                            int ipar = mptr->eval_indices()[index];

                            // Get reference to model parameter
                            const GModelPar& par = (*mptr)[ipar];

                            // Set gradient
                            if (par.is_free()) {
                                if (has_gradient(*mptr, par)) {
                                    (*gradients)[igrad+ipar] =
                                        par.factor_gradient();
                                }
                                else {
                                    (*gradients)[igrad+ipar] =
                                        model_grad(*mptr, par, event);
                                }
                            }
                            else {
                                (*gradients)[igrad+ipar] = 0.0;
                            }

                        } // endfor: looped over all parameter indices

                    } // endif: parameter indices were available

                    // ... otherwise loop over all parameter indices
                    else {

                        // Loop over all parameters
                        for (int ipar = 0; ipar < mptr->size(); ++ipar) {

                            // Get reference to model parameter
                            const GModelPar& par = (*mptr)[ipar];

                            // Set gradient
                            if (par.is_free()) {
                                if (has_gradient(*mptr, par)) {
                                    (*gradients)[igrad+ipar] =
                                        par.factor_gradient();
                                }
                                else {
                                    (*gradients)[igrad+ipar] =
                                        model_grad(*mptr, par, event);
                                }
                            }
                            else {
                                (*gradients)[igrad+ipar] = 0.0;
                            }

                        } // endfor: looped over all parameter indices

                    } // endelse: no parameter indices were available

                } // endif: use model gradients

            } // endif: model component was valid for instrument

            // Increment parameter counter for gradients
            igrad += mptr->size();

        } // endif: model was valid

    } // endfor: Looped over models

    // Return
    return model;
}


/***********************************************************************//**
 * @brief Return vector of model values and (optionally) gradients
 *
 * @param[in] models Model container.
 * @param[out] gradients Pointer to sparse gradient matrix.
 * @return Vector of model values.
 *
 * @exception GException::invalid_argument
 *            Gradient matrix mismatches number of events or model parameters.
 *
 * Returns the model values for each event in the observation.
 *
 * If @p gradients is not NULL, the matrix contains on output the model
 * factor gradients for all events. Each row of the @p gradients matrix
 * corresponds to one event, the columns correspond to the parameters
 * of the @p models container.
 ***************************************************************************/
GVector GObservation::model(const GModels& models,
                            GMatrixSparse* gradients) const
{
    // Initialise variables
    int  nevents  = events()->size();
    bool use_grad = (gradients != NULL);
    int  igrad    = 0;

    // Initialise model values
    GVector values(nevents);

    // If gradient is available then check gradient size and initialise
    // sparse matrix stack
    if (use_grad) {

        // Initialise stack of parameters with gradients
        m_pars_with_gradients.clear();
        m_pars_with_gradients.reserve(models.size());

        // Check number of columns
        int ncolumns = (nevents > 0) ? models.npars() : 0;
        if (gradients->columns() != ncolumns) {
            std::string msg = "Number of "+gammalib::str(gradients->columns())+" "
                              "columns in gradient matrix differs from expected "
                              "number of "+gammalib::str(ncolumns)+". Please "
                              "specify compatible arguments.";
            throw GException::invalid_argument(G_MODEL2, msg);
        }

        // Check number of rows
        int nrows = (ncolumns > 0) ? nevents : 0;
        if (gradients->rows() != nrows) {
            std::string msg = "Number of "+gammalib::str(gradients->rows())+" "
                              "rows in gradient matrix differs from expected "
                              "number of "+gammalib::str(nrows)+". Please "
                              "specify compatible arguments.";
            throw GException::invalid_argument(G_MODEL2, msg);
        }

    } // endif: gradient was requested

    // Loop over models
    for (int i = 0; i < models.size(); ++i) {

        // Get model pointer. Continue only if pointer is valid
        const GModel* mptr = models[i];
        if (mptr != NULL) {

            // Continue only if model applies to specific instrument and
            // observation identifier
            if (mptr->is_valid(instrument(), id())) {

                // If no gradients should be used then evaluate model and
                // add it to values
                if (!use_grad) {
                    values += mptr->eval(*this, NULL);
                }

                // ... otherwise if there are events then evaluate model and
                // gradients and add to model vector and gradient matrix
                else if (nevents > 0) {

                    // Allocate gradient matrix
                    GMatrixSparse grads(nevents, mptr->size());

                    // Initialise sparse matrix stack
                    grads.stack_init(nevents*mptr->size(), mptr->size());

                    // Evaluate model and add to values
                    values += mptr->eval(*this, &grads);

                    // Destroy sparse matrix stack (fills stack in matrix)
                    grads.stack_destroy();

                    // If model provides the parameter indices that were
                    // updated by the eval() method then use them as this
                    // is less time consuming than looping over all indices
                    if (mptr->has_eval_indices()) {

                        // Get number of relevant parameters
                        int npars = mptr->eval_indices().size();

                        // Loop over all relevant parameters
                        for (int index = 0; index < npars; ++index) {

                            // Retrieve parameter index
                            int ipar = mptr->eval_indices()[index];

                            // Get reference to model parameter
                            const GModelPar& par = (*mptr)[ipar];

                            // If parameters is free then set gradient
                            if (par.is_free()) {

                                // Determine gradient
                                GVector grad(nevents);
                                if (has_gradient(*mptr, par)) {
                                    grad = grads.column(ipar);
                                }
                                else {
                                    grad = model_grad(*mptr, par);
                                }

                                // Set gradient
                                gradients->column(igrad+ipar, grad);

                            } // endif: parameter was free

                        } // endfor: looped over all parameter indices

                    } // endif: parameter indices were available

                    // ... otherwise loop over all parameter indices
                    else {

                        // Loop over all parameters
                        for (int ipar = 0; ipar < mptr->size(); ++ipar) {

                            // Get reference to model parameter
                            const GModelPar& par = (*mptr)[ipar];

                            // If parameters is free then set gradient
                            if (par.is_free()) {

                                // Determine gradient
                                GVector grad(nevents);
                                if (has_gradient(*mptr, par)) {
                                    grad = grads.column(ipar);
                                    #if defined(G_DEBUG_VECTOR_MODEL)
                                    if (ipar < 3) {
                                        GVector num_grad = model_grad(*mptr, par);
                                        for (int k = 0; k < nevents; ++k) {
                                            if (grad[k] != 0.0 || num_grad[k] != 0.0) {
                                                std::cout << "ipar=" << ipar;
                                                std::cout << " par=" << par.name();
                                                std::cout << " k=" << k;
                                                std::cout << " ana=" << grad[k];
                                                std::cout << " num=" << num_grad[k];
                                                std::cout << std::endl;
                                            }
                                        }
                                    }
                                    #endif
                                }
                                else {
                                    grad = model_grad(*mptr, par);
                                }

                                // Set gradient
                                gradients->column(igrad+ipar, grad);

                            } // endif: parameter was free

                        } // endfor: looped over all parameter indices

                    } // endelse: no parameter indices were available

                } // endif: use model gradients

            } // endif: model component was valid for instrument

            // Increment parameter counter for gradients
            igrad += mptr->size();

        } // endif: model was valid

    } // endfor: Looped over models

    // Return values
    return values;
}


/***********************************************************************//**
 * @brief Return total number of observed events
 *
 * @returns Total number of observed events.
 *
 * Returns the total number of observed events. If the observation does not
 * contain any events the method returns zero.
 ***************************************************************************/
int GObservation::nobserved(void) const
{
    // Initialise number of observed events
    int nobserved = 0;

    // Extract number of observed events
    if (m_events != NULL) {
        nobserved = m_events->number();
    }

    // Return number of observed events
    return nobserved;
}


/***********************************************************************//**
 * @brief Return total number (and optionally gradients) of predicted counts
 *        for all models
 *
 * @param[in] models Models.
 * @param[out] gradients Model parameter gradients (optional).
 *
 * @exception GException::invalid_argument
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
double GObservation::npred(const GModels& models, GVector* gradients) const
{
    // Verify that gradient vector and models have the same dimension
    #if defined(G_RANGE_CHECK)
    if (gradients != NULL) {
        if (models.npars() != gradients->size()) {
            std::string msg = "Model parameter gradients vector size "+
                              gammalib::str(gradients->size())+" differs from "
                              "number of "+gammalib::str(models.npars())+
                              " model parameters. Please specify a gradients "
                              "vector with the appropriate size.";
            throw GException::invalid_argument(G_NPRED, msg);
        }
    }
    #endif

    // Initialise
    double npred = 0.0;    // Reset predicted number of counts
    int    igrad = 0;      // Reset gradient counter

    // If gradient is available then reset gradient vector elements to 0
    if (gradients != NULL) {
        (*gradients) = 0.0;
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
                npred += this->npred(*mptr);

                // Optionally determine Npred gradients
                if (gradients != NULL) {
                    for (int k = 0; k < mptr->size(); ++k) {
                        const GModelPar& par = (*mptr)[k];
                        (*gradients)[igrad+k] = npred_grad(*mptr, par);
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
 * @brief Return total number of predicted counts for one model
 *
 * @param[in] model Gamma-ray source model.
 *
 * @exception GException::gti_invalid
 *            Good Time Interval is invalid.
 *
 * Computes
 *
 * \f[
 *    N_{\rm pred} = \int_{\rm GTI} \int_{E_{\rm bounds}} \int_{\rm ROI}
 *                   P(p',E',t') \, dp' \, dE' \, dt'
 * \f]
 *
 * of the event probability
 *
 * \f[
 *    P(p',E',t') = \int \int \int
 *                  S(p,E,t) \times R(p',E',t'|p,E,t) \, dp \, dE \, dt
 * \f]
 *
 * where
 * \f$S(p,E,t)\f$ is the source model,
 * \f$R(p',E',t'|p,E,t)\f$ is the instrument response function,
 * \f$p'\f$ is the measured photon direction,
 * \f$E'\f$ is the measured photon energy,
 * \f$t'\f$ is the measured photon arrival time,
 * \f$p\f$ is the true photon arrival direction,
 * \f$E\f$ is the true photon energy, and
 * \f$t\f$ is the true photon arrival time.
 *
 * This method performs the time integration over the Good Time Intervals
 * \f${\rm GTI}\f$. Note that MET is used for the time integration interval.
 * This, however, is no specialisation since npred_grad_kern_spec::eval()
 * converts the argument back in a GTime object by assuming that the argument
 * is in MET, hence the correct time system will be used at the end by the
 * method.
 ***************************************************************************/
double GObservation::npred(const GModel& model) const
{
    // Initialise result
    double npred = 0.0;

    // Continue only if model applies to specific instrument and
    // observation identifier
    if (model.is_valid(instrument(), id())) {

        // Case A: If the model is constant then integrate analytically
        if (model.is_constant()) {

            // Evaluate model at first start time and multiply by ontime
            double ontime = events()->gti().ontime();

            // Integrate only if ontime is positive
            if (ontime > 0.0) {

                // Integration is a simple multiplication by the time
                npred = npred_spec(model, events()->gti().tstart()) * ontime;

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
                    std::string msg = "Start time "+
                                      events()->gti().tstart(i).print()+
                                      " is past the stop time "+
                                      events()->gti().tstop(i).print()+
                                      " for Good Time Interval "+
                                      gammalib::str(i)+". Please make sure "
                                      "that all Good Time Intervals have start "
                                      "times that are earlier than the stop "
                                      "time.";
                    throw GException::invalid_value(G_NPRED, msg);
                }

                // Setup integration function
                GObservation::npred_kern integrand(this, &model);
                GIntegral                integral(&integrand);

                // Do Romberg integration
                npred += integral.romberg(tstart, tstop);

            } // endfor: looped over GTIs

        } // endelse: integrated temporally

    } // endif: model was valid

    // Return Npred
    return npred;
}


/***********************************************************************//**
 * @brief Returns parameter gradient of model for a given event
 *
 * @param[in] model Model.
 * @param[in] par Model parameter.
 * @param[in] event Event.
 *
 * This method uses a robust but simple difference method to estimate
 * parameter gradients that have not been provided by the model. We use here
 * a simple method as this method is likely used for spatial model parameters,
 * and the spatial model may eventually be noisy due to numerical integration
 * limits.
 *
 * The step size for the simple method has been fixed to 0.0002, which
 * corresponds to about 1 arcsec for parameters that are given in degrees.
 * The reasoning behind this value is that parameters that use numerical
 * gradients are typically angles, such as for example the position, and
 * we want to achieve arcsec precision with this method.
 ***************************************************************************/
double GObservation::model_grad(const GModel&    model,
                                const GModelPar& par,
                                const GEvent&    event) const
{
    // Initialise gradient
    double grad = 0.0;

    // Compute gradient only if parameter is free
    if (par.is_free()) {

        // Get non-const model pointer
        GModelPar* ptr = const_cast<GModelPar*>(&par);

        // Save current model parameter
        GModelPar current = par;

        // Get actual parameter value
        double x = par.factor_value();

        // Set fixed step size for computation of derivative.
        // By default, the step size is fixed to 0.0002.
        const double step_size = 0.0002; // ~1 arcsec
        double       h         = step_size;

        // Re-adjust the step-size h in case that the initial step size is
        // larger than the allowed parameter range 
        if (par.has_min() && par.has_max()) {
            double par_h = par.factor_max() - par.factor_min();
            if (par_h < h) {
                h = par_h;
            }
        }

        // Continue only if step size is positive
        if (h > 0.0) {

            // Remove any boundaries to avoid limitations
            //ptr->remove_range(); // Not needed in principle

            // Setup derivative function
            GObservation::model_func function(this, &model, ptr, &event);
            GDerivative              derivative(&function);

            // If we are too close to the minimum boundary use a right sided
            // difference ...
            if (par.has_min() && ((x-par.factor_min()) < h)) {
                grad = derivative.right_difference(x, h);
            }

            // ... otherwise if we are too close to the maximum boundary use
            // a left sided difference ...
            else if (par.has_max() && ((par.factor_max()-x) < h)) {
                grad = derivative.left_difference(x, h);
            }

            // ... otherwise use a symmetric difference
            else {
                grad = derivative.difference(x, h);
            }

        } // endif: step size was positive

        // Restore current model parameter
        *ptr = current;

    } // endif: model parameter was free

    // Return gradient
    return grad;
}


/***********************************************************************//**
 * @brief Returns parameter gradients of model for all events
 *
 * @param[in] model Model.
 * @param[in] par Model parameter.
 *
 * This method uses a robust but simple difference method to estimate
 * parameter gradients that have not been provided by the model. We use here
 * a simple method as this method is likely used for spatial model parameters,
 * and the spatial model may eventually be noisy due to numerical integration
 * limits.
 *
 * The step size for the simple method has been fixed to 0.0002, which
 * corresponds to about 1 arcsec for parameters that are given in degrees.
 * The reasoning behind this value is that parameters that use numerical
 * gradients are typically angles, such as for example the position, and
 * we want to achieve arcsec precision with this method.
 ***************************************************************************/
GVector GObservation::model_grad(const GModel&    model,
                                 const GModelPar& par) const
{
    // Initialise gradients
    GVector gradients(events()->size());

    // Compute gradient only if parameter is free
    if (par.is_free()) {

        // Get non-const model pointer
        GModelPar* ptr = const_cast<GModelPar*>(&par);

        // Save current model parameter
        GModelPar current = par;

        // Get actual parameter value
        double x = par.factor_value();

        // Set fixed step size for computation of derivative.
        // By default, the step size is fixed to 0.0002.
        const double step_size = 0.0002; // ~1 arcsec
        double       h         = step_size;

        // Re-adjust the step-size h in case that the initial step size is
        // larger than the allowed parameter range
        if (par.has_min() && par.has_max()) {
            double par_h = par.factor_max() - par.factor_min();
            if (par_h < h) {
                h = par_h;
            }
        }

        // Continue only if step size is positive
        if (h > 0.0) {

            // Setup derivative function
            //GObservation::model_func function(this, &model, ptr, &event);
            //GDerivative              derivative(&function);

            // If we are too close to the minimum boundary use a right sided
            // difference ...
            if (par.has_min() && ((x-par.factor_min()) < h)) {
                //grad = derivative.right_difference(x, h);
                // Compute fs1
                ptr->factor_value(x);
                GVector fs1 = model.eval(*this, NULL);

                // Compute fs2
                ptr->factor_value(x+h);
                GVector fs2 = model.eval(*this, NULL);

                // Compute derivative
                gradients = (fs1 - fs2) / h;

            }

            // ... otherwise if we are too close to the maximum boundary use
            // a left sided difference ...
            else if (par.has_max() && ((par.factor_max()-x) < h)) {
                //grad = derivative.left_difference(x, h);
                // Compute fs1
                ptr->factor_value(x);
                GVector fs1 = model.eval(*this, NULL);

                // Compute fs2
                ptr->factor_value(x-h);
                GVector fs2 = model.eval(*this, NULL);

                // Compute derivative
                gradients = (fs1 - fs2) / h;
            }

            // ... otherwise use a symmetric difference
            else {
                //grad = derivative.difference(x, h);
                // Compute fs1
                ptr->factor_value(x+h);
                GVector fs1 = model.eval(*this, NULL);

                // Compute fs2
                ptr->factor_value(x-h);
                GVector fs2 = model.eval(*this, NULL);

                // Compute derivative
                gradients = 0.5 * (fs1 - fs2) / h;
            }

        } // endif: step size was positive

        // Restore current model parameter
        *ptr = current;

    } // endif: model parameter was free

    // Return gradients
    return gradients;
}


/***********************************************************************//**
 * @brief Returns parameter gradient of Npred
 *
 * @param[in] model Gamma-ray source model.
 * @param[in] par Model parameter.
 *
 * Computes
 *
 * \f[
 *    \frac{{\rm d} N_{\rm pred}}{{\rm d} a_i}
 * \f]
 *
 * where
 *
 * \f[
 *    N_{\rm pred} = \int_{\rm GTI} \int_{E_{\rm bounds}} \int_{\rm ROI}
 *                   P(p',E',t') \, dp' \, dE' \, dt'
 * \f]
 *
 * is the integral of the event probability
 *
 * \f[
 *    P(p',E',t') = \int \int \int
 *                  S(p,E,t) \times R(p',E',t'|p,E,t) \, dp \, dE \, dt
 * \f]
 *
 * and
 * \f$a_i\f$ is the model parameter \f$i\f$.
 * Furthermore
 * \f$S(p,E,t)\f$ is the source model,
 * \f$R(p',E',t'|p,E,t)\f$ is the instrument response function,
 * \f$p'\f$ is the measured photon direction,
 * \f$E'\f$ is the measured photon energy,
 * \f$t'\f$ is the measured photon arrival time,
 * \f$p\f$ is the true photon arrival direction,
 * \f$E\f$ is the true photon energy, and
 * \f$t\f$ is the true photon arrival time.
 *
 * This method uses a robust but simple difference method to estimate
 * parameter gradients that have not been provided by the model. We use here
 * a simple method as this method is likely used for spatial model parameters,
 * and the spatial model may eventually be noisy due to numerical integration
 * limits.
 *
 * The step size for the simple method has been fixed to 0.0002.
 ***************************************************************************/
double GObservation::npred_grad(const GModel& model, const GModelPar& par) const
{
    // Initialise result
    double grad = 0.0;

    // Compute gradient only if parameter is free
    if (par.is_free()) {

        // Get non-const model pointer
        GModelPar* ptr = const_cast<GModelPar*>(&par);

        // Save current model parameter
        GModelPar current = par;

        // Get actual parameter value
        double x = par.factor_value();

        // Set fixed step size to 0.0002 for computation of derivative.
        const double step_size = 0.0002;
        double       h         = step_size;

        // Re-adjust the step-size h in case that the initial step size is
        // larger than the allowed parameter range 
        if (par.has_min() && par.has_max()) {
            double par_h = par.factor_max() - par.factor_min();
            if (par_h < h) {
                h = par_h;
            }
        }

        // Continue only if step size is positive
        if (h > 0.0) {

            // Remove any boundaries to avoid limitations
            //ptr->remove_range(); // Not needed in principle

            // Setup derivative function
            GObservation::npred_func function(this, &model, ptr);
            GDerivative              derivative(&function);

            // If we are too close to the minimum boundary use a right sided
            // difference ...
            if (par.has_min() && ((x-par.factor_min()) < h)) {
                grad = derivative.right_difference(x, h);
            }

            // ... otherwise if we are too close to the maximum boundary use
            // a left sided difference ...
            else if (par.has_max() && ((par.factor_max()-x) < h)) {
                grad = derivative.left_difference(x, h);
            }

            // ... otherwise use a symmetric difference
            else {
                grad = derivative.difference(x, h);
            }

        } // endif: step size was positive

        // Restore current model parameter
        *ptr = current;

    } // endif: model parameter was free

    // Return result
    return grad;
}


/***********************************************************************//**
 * @brief Response cache removal hook
 *
 * @param[in] name Model name.
 *
 * Remove response cache for model @p name from response cache.
 ***************************************************************************/
void GObservation::remove_response_cache(const std::string& name)
{
    // Build model name
    std::string model_name  = id() + ":" + name;

    // Remove response cache
    const_cast<GResponse*>(this->response())->remove_response_cache(model_name);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Check whether a model parameter has an analytical gradient
 *
 * @param[in] model Model.
 * @param[in] par Model parameter.
 * @returns True if model parameter is free and has an analytical gradient.
 *
 * Checks whether a model parameter is free and has an analytical gradient.
 ***************************************************************************/
bool GObservation::has_gradient(const GModel& model, const GModelPar& par) const
{
    // Initialise flag
    bool has_gradient = false;

    // Only consider free parameters
    if (par.is_free()) {

        // Build identifier
        std::string id = model.name() + ":" + par.name();

        // Search for parameter address in stack and set flag to true if
        // parameter address was found
        for (int i = 0; i < m_pars_with_gradients.size(); ++i) {
            if (m_pars_with_gradients[i] == id) {
                has_gradient = true;
                break;
            }
        }

    } // endif: parameter was free

    // Return flag
    return has_gradient;
}


/***********************************************************************//**
 * @brief Signals that an analytical gradient was computed for a model
 *        parameter
 *
 * @param[in] model Model.
 * @param[in] par Model parameter.
 *
 * Signals that an analytical gradient was computed for a model parameter.
 ***************************************************************************/
void GObservation::computed_gradient(const GModel& model, const GModelPar& par) const
{
    // Initialise flag
    bool in_stack = false;

    // Build identifier
    std::string id = model.name() + ":" + par.name();

    // Check if parameter is already in stack
    for (int i = 0; i < m_pars_with_gradients.size(); ++i) {
        if (m_pars_with_gradients[i] == id) {
            in_stack = true;
            break;
        }
    }

    // If parameter is not yet in stack the push it on the stack
    if (!in_stack) {
        m_pars_with_gradients.push_back(id);
    }

    // Return
    return;
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
    m_statistic = "cstat";
    m_events    = NULL;

    // Initialise stack of identifiers of parameters with gradients
    m_pars_with_gradients.clear();

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
    m_name      = obs.m_name;
    m_id        = obs.m_id;
    m_statistic = obs.m_statistic;

    // Copy stack of identifiers of parameters with gradients
    m_pars_with_gradients = obs.m_pars_with_gradients;

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
 * @brief Evaluate log-likelihood function for Poisson statistic and
 *        unbinned analysis (version with working arrays)
 *
 * @param[in] models Models.
 * @param[in,out] gradients Gradient.
 * @param[in,out] curvature Curvature matrix.
 * @param[in,out] npred Number of predicted events.
 * @return Likelihood value.
 *
 * This method evaluates the -(log-likelihood) function for parameter
 * optimisation using unbinned analysis and Poisson statistic.
 * The -(log-likelihood) function is given by
 *
 * \f[
 *    L = N_{\rm pred} - \sum_i \log e_i
 * \f]
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
                                                 GVector*       gradients,
                                                 GMatrixSparse* curvature,
                                                 double*        npred) const
{
    // Initialise likelihood value
    double value = 0.0;

    // Get number of events and parameters
    int nevents = events()->size();
    int npars   = gradients->size();

    // Allocate some working arrays
    int*          inx    = new int[npars];
    double*       values = new double[npars];
    GMatrixSparse wrk_matrix(nevents, npars);
    GVector       wrk_grad(npars);

    // Determine Npred value and gradient for this observation
    double npred_value = this->npred(models, &wrk_grad);

    // Update likelihood, Npred and gradients
    value      += npred_value;
    *npred     += npred_value;
    *gradients += wrk_grad;

    // Compute model and derivative
    GVector model_vector = this->model(models, &wrk_matrix);

    // Transpose matrix
    wrk_matrix = wrk_matrix.transpose();

    // Iterate over all events
    for (int i = 0; i < nevents; ++i) {

        // Get event pointer
        const GEvent* event = (*events())[i];

        // Get model value
        double model = model_vector[i];

        // Skip bin if model is too small (avoids -Inf or NaN gradients)
        if (model <= minmod) {
            continue;
        }

        // Extract working gradient multiplied by bin size
        GVector wrk_grad = wrk_matrix.column(i);

        // Create index array of non-zero derivatives and initialise working
        // array
        int ndev = 0;
        for (int ipar = 0; ipar < npars; ++ipar) {
            values[ipar] = 0.0;
            if (wrk_grad[ipar] != 0.0 && !gammalib::is_infinite(wrk_grad[ipar])) {
                inx[ndev] = ipar;
                ndev++;
            }
        }

        // Update Poissonian statistic (excluding factorial term for faster
        // computation)
        value -= std::log(model);

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
            (*gradients)[jpar] -= fb * g;

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
 * @brief Evaluate log-likelihood function for Poisson statistic and
 *        binned analysis (version with working arrays)
 *
 * @param[in] models Models.
 * @param[in,out] gradients Gradient.
 * @param[in,out] curvature Curvature matrix.
 * @param[in,out] npred Number of predicted events.
 * @return Likelihood value.
 *
 * This method evaluates the -(log-likelihood) function for parameter
 * optimisation using binned analysis and Poisson statistic.
 * The -(log-likelihood) function is given by
 *
 * \f[
 *    L=-\sum_i n_i \log e_i - e_i
 * \f]
 *
 * where the sum is taken over all data space bins, \f$n_i\f$ is the
 * observed number of counts and \f$e_i\f$ is the model.
 * This method also computes the parameter gradients
 * \f$\delta L/dp\f$
 * and the curvature matrix
 * \f$\delta^2 L/dp_1 dp_2\f$
 * and also updates the total number of predicted events m_npred.
 ***************************************************************************/
double GObservation::likelihood_poisson_binned(const GModels& models,
                                               GVector*       gradients,
                                               GMatrixSparse* curvature,
                                               double*        npred) const
{
    // Initialise likelihood value
    double value = 0.0;

    // Initialise statistic
    #if defined(G_OPT_DEBUG)
    int    n_bins         = 0;
    int    n_used         = 0;
    int    n_small_model  = 0;
    int    n_zero_data    = 0;
    int    n_exclude_data = 0;
    double sum_data       = 0.0;
    double sum_model      = 0.0;
    double init_value     = value;
    #endif

    // Get number of events and parameters
    int nevents = events()->size();
    int npars   = gradients->size();

    // Allocate some working arrays
    int*          inx    = new int[npars];
    double*       values = new double[npars];
    GMatrixSparse wrk_matrix(nevents, npars);

    // Compute model and derivative
    GVector model_vector = this->model(models, &wrk_matrix);

    // Transpose matrix
    wrk_matrix = wrk_matrix.transpose();

    // Iterate over all bins
    for (int i = 0; i < nevents; ++i) {

        // Update number of bins
        #if defined(G_OPT_DEBUG)
        n_bins++;
        #endif

        // Get event pointer
        const GEventBin* bin =
            (*(static_cast<GEventCube*>(const_cast<GEvents*>(events()))))[i];

        // Get number of counts in bin
        double data = bin->counts();

        // Skip bin if data is negative (filtering flag)
        if (data < 0) {
            #if defined(G_OPT_DEBUG)
            n_exclude_data++;
            #endif
            continue;
        }

        // Get model value
        double model = model_vector[i] * bin->size();

        // Skip bin if model is too small (avoids -Inf or NaN gradients)
        if (model <= minmod) {
            #if defined(G_OPT_DEBUG)
            n_small_model++;
            #endif
            continue;
        }

        // Update statistic
        #if defined(G_OPT_DEBUG)
        n_used++;
        sum_data  += data;
        sum_model += model;
        #endif

        // Update Npred
        *npred += model;

        // Extract working gradient multiplied by bin size
        GVector wrk_grad = wrk_matrix.column(i) * bin->size();

        // Create index array of non-zero derivatives and initialise working
        // array
        int ndev = 0;
        for (int ipar = 0; ipar < npars; ++ipar) {
            values[ipar] = 0.0;
            if (wrk_grad[ipar] != 0.0 && !gammalib::is_infinite(wrk_grad[ipar])) {
                inx[ndev] = ipar;
                ndev++;
            }
        }

        // Update gradient vector and curvature matrix. To avoid
        // unneccessary computations we distinguish the case where
        // data>0 and data=0. The second case requires much less
        // computation since it does not contribute to the curvature
        // matrix ...
        if (data > 0.0) {

            // Update Poissonian statistic (excluding factorial term for
            // faster computation)
            value -= data * std::log(model) - model;

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
                (*gradients)[jpar] += fc * g;

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

            // Update statistic
            #if defined(G_OPT_DEBUG)
            n_zero_data++;
            #endif

            // Update Poissonian statistic (excluding factorial term for
            // faster computation)
            value += model;

            // Skip bin now if there are no non-zero derivatives
            if (ndev < 1) {
                continue;
            }

            // Update gradient
            register int* ipar = inx;
            for (register int idev = 0; idev < ndev; ++idev, ++ipar) {
                (*gradients)[*ipar] += wrk_grad[*ipar];
            }

        } // endif: data was 0

    } // endfor: iterated over all events

    // Free temporary memory
    if (values != NULL) delete [] values;
    if (inx    != NULL) delete [] inx;

    // Dump statistic
    #if defined(G_OPT_DEBUG)
    std::cout << "Number of bins: " << n_bins << std::endl;
    std::cout << "Number of bins used for computation: " << n_used << std::endl;
    std::cout << "Number of bins excluded from data: " << n_exclude_data << std::endl;
    std::cout << "Number of bins excluded due to small model: " << n_small_model << std::endl;
    std::cout << "Number of bins with zero data: " << n_zero_data << std::endl;
    std::cout << "Sum of data: " << sum_data << std::endl;
    std::cout << "Sum of model: " << sum_model << std::endl;
    std::cout << "Initial statistic: " << init_value << std::endl;
    std::cout << "Statistic: " << value-init_value << std::endl;
    #endif

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate log-likelihood function for Gaussian statistic and
 *        binned analysis (version with working arrays)
 *
 * @param[in] models Models.
 * @param[in,out] gradients Gradient.
 * @param[in,out] curvature Curvature matrix.
 * @param[in,out] npred Number of predicted events.
 * @return Likelihood value.
 *
 * This method evaluates the -(log-likelihood) function for parameter
 * optimisation using binned analysis and Poisson statistic.
 * The -(log-likelihood) function is given by
 *
 * \f[
 *    L = 1/2 \sum_i (n_i - e_i)^2 \sigma_i^{-2}
 * \f]
 *
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
                                                GVector*       gradients,
                                                GMatrixSparse* curvature,
                                                double*        npred) const
{
    // Initialise likelihood value
    double value = 0.0;

    // Get number of events and parameters
    int nevents = events()->size();
    int npars   = gradients->size();

    // Allocate some working arrays
    int*          inx    = new int[npars];
    double*       values = new double[npars];
    GMatrixSparse wrk_matrix(nevents, npars);

    // Compute model and derivative
    GVector model_vector = this->model(models, &wrk_matrix);

    // Transpose matrix
    wrk_matrix = wrk_matrix.transpose();

    // Iterate over all bins
    for (int i = 0; i < nevents; ++i) {

        // Get event pointer
        const GEventBin* bin =
            (*(static_cast<GEventCube*>(const_cast<GEvents*>(events()))))[i];

        // Get number of counts in bin
        double data = bin->counts();

        // Skip bin if data is negative (filtering flag)
        if (data < 0) {
            continue;
        }

        // Get statistical uncertainty
        double sigma = bin->error();

        // Skip bin if statistical uncertainty is too small
        if (sigma <= minerr) {
            continue;
        }

        // Get model value
        double model = model_vector[i] * bin->size();

        // Skip bin if model is too small (avoids -Inf or NaN gradients)
        if (model <= minmod) {
            continue;
        }

        // Update Npred
        *npred += model;

        // Extract working gradient multiplied by bin size
        GVector wrk_grad = wrk_matrix.column(i) * bin->size();

        // Create index array of non-zero derivatives and initialise working
        // array
        int ndev = 0;
        for (int ipar = 0; ipar < npars; ++ipar) {
            values[ipar] = 0.0;
            if (wrk_grad[ipar] != 0.0 && !gammalib::is_infinite(wrk_grad[ipar])) {
                inx[ndev] = ipar;
                ndev++;
            }
        }

        // Set weight
        double weight = 1.0 / (sigma * sigma);

        // Update Gaussian statistic
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
            (*gradients)[jpar] -= fa * fa_i;

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
    // Set value
    m_par->factor_value(x);

    // Compute model value
    double value = m_model->eval(*m_event, *m_parent);

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
    // Set value
    m_par->factor_value(x);

    // Compute Npred value
    double npred = m_parent->npred(*m_model);

    // Return value
    return npred;
}


/***********************************************************************//**
 * @brief Integration kernel for npred() method
 *
 * @param[in] x Function value.
 ***************************************************************************/
double GObservation::npred_kern::eval(const double& x)
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
 *
 * \f[
 *    N_{\rm pred}(t') = \int_{E_{\rm bounds}} \int_{\rm ROI}
 *                       P(p',E',t') \, dp' \, dE'
 * \f]
 *
 * of the event probability
 *
 * \f[
 *    P(p',E',t') = \int \int \int
 *                  S(p,E,t) \times R(p',E',t'|p,E,t) \, dp \, dE \, dt
 * \f]
 *
 * where
 * \f$S(p,E,t)\f$ is the source model,
 * \f$R(p',E',t'|p,E,t)\f$ is the instrument response function,
 * \f$p'\f$ is the measured photon direction,
 * \f$E'\f$ is the measured photon energy,
 * \f$t'\f$ is the measured photon arrival time,
 * \f$p\f$ is the true photon arrival direction,
 * \f$E\f$ is the true photon energy, and
 * \f$t\f$ is the true photon arrival time.
 *
 * This method performs the energy intergration over the energy boundaries
 * \f$E_{\rm bounds}\f$.
 ***************************************************************************/
double GObservation::npred_spec(const GModel& model,
                                const GTime&  obsTime) const
{
    // Set number of iterations for Romberg integration.
    static const int iter = 8;

    // Initialise result
    double result = 0.0;

    // Get energy boundaries
    GEbounds ebounds = events()->ebounds();

    // Loop over energy boundaries
    for (int i = 0; i < ebounds.size(); ++i) {

        // Get boundaries in MeV
        double emin = ebounds.emin(i).MeV();
        double emax = ebounds.emax(i).MeV();

        // Continue only if valid
        if (emax > emin) {

            // Setup integration function
            GObservation::npred_spec_kern integrand(this, &model, &obsTime);
            GIntegral                     integral(&integrand);

            // Set number of iterations
            integral.fixed_iter(iter);

            // Do Romberg integration
            #if defined(G_LN_ENERGY_INT)
            emin = std::log(emin);
            emax = std::log(emax);
            #endif
            result += integral.romberg(emin, emax);

        } // endif: energy interval was valid

    } // endfor: looped over energy boundaries

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(result) || gammalib::is_infinite(result)) {
        std::cout << "*** ERROR: GObservation::npred_spec:";
        std::cout << " NaN/Inf encountered";
        std::cout << " (result=" << result;
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
