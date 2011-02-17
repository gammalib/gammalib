/***************************************************************************
 *           GObservation.cpp  -  Abstract observation base class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GObservation.cpp
 * @brief Abstract observation base class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GObservation.hpp"
#include "GModelSky.hpp"
#include "GModelData.hpp"
#include "GModelSpatialPtsrc.hpp"
#include "GIntegral.hpp"
#include "GDerivative.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
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

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_LN_ENERGY_INT  1   //!< ln(E) variable substitution for integration

/* __ Debug definitions __________________________________________________ */


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
 *
 * Assign observation.
 ***************************************************************************/
GObservation& GObservation::operator= (const GObservation& obs)
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
 * parameter derivative of this probability.
 *
 * If NULL is passed for the gradient vector, then gradients will not be
 * computed.
 ***************************************************************************/
double GObservation::model(const GModels& models, const GEvent& event,
                           GVector* gradient) const
{
    // Verify that gradient vector and models have the same dimension
    #if defined(G_RANGE_CHECK)
    if (gradient != NULL) {
        if (models.npars() != gradient->size())
            throw GException::gradient_par_mismatch(G_MODEL, 
                                                    gradient->size(),
                                                    models.npars());
    }
    #endif

    // Initialise
    double model = 0.0;    // Reset model value
    int    igrad = 0;      // Reset gradient counter

    // If gradient is available then reset gradient vector elements to 0
    if (gradient != NULL)
        (*gradient) = 0.0;

    // Loop over models
    for (int i = 0; i < models.size(); ++i) {

        // Check if model applies to specific instrument
        if (models[i].isvalid(instrument())) {

            // Compute value and add to model
            model += models[i].eval_gradients(event, *this);

            // Optionally determine model gradients
            if (gradient != NULL) {
                for (int k = 0; k < models[i].size(); ++k)
                    (*gradient)[igrad+k] = model_grad(models[i], event, k);
            }

        } // endif: model component was valid for instrument

        // Increment parameter counter for gradients
        igrad += models[i].size();

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
 ***************************************************************************/
double GObservation::npred(const GModels& models, GVector* gradient) const
{
    // Verify that gradient vector and models have the same dimension
    #if defined(G_RANGE_CHECK)
    if (gradient != NULL) {
        if (models.npars() != gradient->size())
            throw GException::gradient_par_mismatch(G_MODEL, 
                                                    gradient->size(),
                                                    models.npars());
    }
    #endif

    // Initialise
    double npred = 0.0;    // Reset predicted number of counts
    int    igrad = 0;      // Reset gradient counter

    // If gradient is available then reset gradient vector elements to 0
    if (gradient != NULL)
        (*gradient) = 0.0;

    // Loop over models
    for (int i = 0; i < models.size(); ++i) {

        // Handle only components that are relevant for the actual
        // instrument
        if (models[i].isvalid(instrument())) {

            // Determine Npred for model
            npred += npred_temp(models[i]);

            // Optionally determine Npred gradients
            if (gradient != NULL) {
                for (int k = 0; k < models[i].size(); ++k)
                    (*gradient)[igrad+k] = npred_grad(models[i], k);
            }

        } // endif: model component was valid for instrument

        // Increment parameter counter for gradient
        igrad += models[i].size();

    } // endfor: Looped over models

    // Return prediction
    return npred;
}


/***********************************************************************//**
 * @brief Set observation name
 *
 * @param[in] name Observation name.
 *
 * Set name of the observation.
 ***************************************************************************/
void GObservation::name(const std::string& name)
{
    // Set name
    m_name = name;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set event container
 *
 * @param[in] events Event container pointer.
 *
 * Set the event container for this observation by cloning the container
 * specified in the argument. If NULL is passed to this method, any existing
 * events are cleared an no event container is attached.
 ***************************************************************************/
void GObservation::events(const GEvents* events)
{
    // Remove an existing event container
    if (m_events != NULL) delete m_events;

    // Signal event container as free
    m_events = NULL;

    // Set event container if the input pointer is valid
    if (events != NULL)
        m_events = events->clone();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set optimizer statistics
 *
 * @param[in] statistics Optimizer statistics.
 *
 * Set optimizer statistics for the observation.
 ***************************************************************************/
void GObservation::statistics(const std::string& statistics)
{
    // Set statistics
    m_statistics = statistics;

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
    if (m_events == NULL)
        throw GException::no_events(G_EVENTS);

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
 =                        Npred integration methods                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Returns parameter gradient of model for a given event
 *
 * @param[in] model Model.
 * @param[in] event Event.
 * @param[in] ipar Parameter index for which gradient should be returned.
 ***************************************************************************/
double GObservation::model_grad(const GModel& model, const GEvent& event,
                                int ipar) const
{
    // Initialise gradient
    double grad = 0.0;

    // Compute gradient only if parameter is free
    if (model[ipar].isfree()) {

        // If model has a gradient then use it
        if (model[ipar].hasgrad())
            grad = model[ipar].gradient();

        // ... otherwise compute it numerically
        else {

            // Setup derivative function
            GObservation::model_func function(this, model, event, ipar);

            // Get derivative
            GDerivative derivative(&function);
            double x = model[ipar].value();
            derivative.eps(1.0e-4);
            derivative.silent(true);
            grad = derivative.value(x, 0.2*std::abs(x));

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
 ***************************************************************************/
double GObservation::npred_grad(const GModel& model, int ipar) const
{
    // Initialise result
    double result = 0.0;

    // Compute gradient only if parameter is free
    if (model[ipar].isfree()) {

        // Setup derivative function
        GObservation::npred_func function(this, model, ipar);
        GDerivative              derivative(&function);

        // Get derivative
        result = derivative.value(model[ipar].value());
        
    } // endif: model parameter was free

    // Return result
    return result;
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

    // Loop over GTIs
    for (int i = 0; i < events()->gti().size(); ++i) {

        // Set integration interval in MET
        double tstart = events()->gti().tstart(i).met();
        double tstop  = events()->gti().tstop(i).met();

        // Throw exception if time interval is not valid
        if (tstop <= tstart)
            throw GException::gti_invalid(G_NPRED_TEMP, events()->gti().tstart(i),
                                          events()->gti().tstop(i));

        // Setup integration function
        GObservation::npred_kern_spec integrand(this, model);
        GIntegral                     integral(&integrand);

        // Do Romberg integration
        result += integral.romb(tstart, tstop);

    } // endfor: looped over GTIs

    // Return result
    return result;
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
    double emin = events()->ebounds().emin().MeV();
    double emax = events()->ebounds().emax().MeV();

    // Throw exception if energy range is not valid
    if (emax <= emin)
        throw GException::erange_invalid(G_NPRED_SPEC, emin, emax);

    // Setup integration function
    GObservation::npred_kern_spat integrand(this, model, obsTime);
    GIntegral                     integral(&integrand);


    // Do Romberg integration
    #if G_LN_ENERGY_INT
    emin = log(emin);
    emax = log(emax);
    #endif
    double result = integral.romb(emin, emax);

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Integration kernel for npred_spec() method
 *
 * @param[in] x Function value.
 *
 * This method implements the integration kernel needed for the npred_spec()
 * method. Upon the defintion of the G_LN_ENERGY_INT declaration the energy
 * integration is done logarithmically (G_LN_ENERGY_INT=1) or not.
 ***************************************************************************/
double GObservation::npred_kern_spat::eval(double x)
{
    #if G_LN_ENERGY_INT
    // Variable substitution
    x = exp(x);
    #endif

    // Set energy in MeV
    GEnergy eng;
    eng.MeV(x);

    // Get function value
    double value = m_model->npred(eng, *m_time, *m_parent);

    #if G_LN_ENERGY_INT
    // Correct for variable substitution
    value *= x;
    #endif

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Integration kernel for npred_temp() method
 *
 * @param[in] x Function value.
 *
 * Note that MET is used for the time conversion. This, however, is no
 * specialisation since npred_grad_temp() hands MET.
 ***************************************************************************/
double GObservation::npred_kern_spec::eval(double x)
{
    // Convert argument in MET
    GTime time;
    time.met(x);

    // Return value
    return (m_parent->npred_spec(*m_model,time));
}


/*==========================================================================
 =                                                                         =
 =                        Model gradient computation                       =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Model function evaluation for gradient computation
 *
 * @param[in] x Function value.
 *
 * @todo We simply remove any parameter boundaries here for the computation
 *       to avoid any out of boundary errors. We may have models, however,
 *       for which out of bound parameters lead to illegal computations, such
 *       a division by zero or taking the square root of negative values.
 *       I cannot see any elegant method to catch this at this level.
 *       Eventually, the higher level method should avoid going in a
 *       parameter domain that is not defined. 
 ***************************************************************************/
double GObservation::model_func::eval(double x)
{
    // Circumvent const correctness
    GModel* model = (GModel*)m_model;
    
    // Save current model parameter
    GModelPar current = (*model)[m_ipar];

    // Set requested model value. Remove any boundaries to avoid limitations.
    (*model)[m_ipar].value(x);
    (*model)[m_ipar].remove_range();

    // Compute model value
    double value = model->eval_gradients(*m_event, *m_parent);

    // Restore current model parameter
    (*model)[m_ipar] = current;

    // Return value
    return value;
}


/*==========================================================================
 =                                                                         =
 =                    Npred gradient integration methods                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Npred function evaluation for gradient computation
 *
 * @param[in] x Function value.
 *
 * @todo We simply remove any parameter boundaries here for the computation
 *       to avoid any out of boundary errors. We may have models, however,
 *       for which out of bound parameters lead to illegal computations, such
 *       a division by zero or taking the square root of negative values.
 *       I cannot see any elegant method to catch this at this level.
 *       Eventually, the higher level method should avoid going in a
 *       parameter domain that is not defined. 
 ***************************************************************************/
double GObservation::npred_func::eval(double x)
{
    // Circumvent const correctness
    GModel* model = (GModel*)m_model;
    
    // Save current model parameter
    GModelPar current = (*model)[m_ipar];

    // Set requested model value. Remove any boundaries to avoid limitations.
    (*model)[m_ipar].value(x);
    (*model)[m_ipar].remove_range();

    // Compute Npred value
    double npred = m_parent->npred_temp(*model);

    // Restore current model parameter
    (*model)[m_ipar] = current;

    // Return value
    return npred;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] obs Observation.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GObservation& obs)
{
     // Write observation in output stream
    os << obs.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] obs Observation.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GObservation& obs)
{
    // Write observation into logger
    log << obs.print();

    // Return logger
    return log;
}
