/***************************************************************************
 *           GObservation.cpp  -  Observation abstract base class          *
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
 * @brief GObservation abstract base class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GException.hpp"
#include "GObservation.hpp"
#include "GModelSky.hpp"
#include "GModelSpatialPtsrc.hpp"
#include "GIntegral.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_MODEL                   "GObservation::model(GModels&, GPointing&,"\
                                    " GInstDir&, GEnergy&, GTime&, GVector*)"
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
#define G_LN_ENERGY_INT 1    //!< ln(E) variable substitution for integration

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GObservation::GObservation(void)
{
    // Initialise class members for clean destruction
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
    // Initialise class members for clean destruction
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
 * Assign one observation to another.
 ***************************************************************************/
GObservation& GObservation::operator= (const GObservation& obs)
{
    // Execute only if object is not identical
    if (this != &obs) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
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
 * @param[out] gradient Pointer to gradient vector (NULL=not computed).
 *
 * @exception GException::gradient_par_mismatch
 *            Dimension of gradient vector mismatches number of parameters.
 *
 * Implements generic model and gradient evaluation for an observation. The
 * model gives the probability for an event to occur with a given instrument
 * direction, at a given energy and at a given time. The gradient is the
 * parameter derivative of this probability.
 *
 * @todo We actually circumvent the const correctness for the eval_gradients()
 * method, yet I don't know why this method is not const. This should be
 * checked, and better put this method const.
 ***************************************************************************/
double GObservation::model(const GModels& models, const GEvent& event,
                           GVector* gradient) const
{
    // Verify that gradients vector has the same dimension than the
    // model has parameters
    #if defined(G_RANGE_CHECK)
    if (models.npars() != gradient->size())
        throw GException::gradient_par_mismatch(G_MODEL, gradient->size(),
                                                models.npars());
    #endif

    // Initialise model value and gradient index
    double model = 0.0;
    int    igrad = 0;

    // Loop over models
    for (int i = 0; i < models.size(); ++i) {

        // Get model pointer (enforce const correctness)
        GModel* ptr = (GModel*)models(i);

        // Check if model applies to specific instrument
        if (ptr->isvalid(instrument())) {

            // Compute value and add to model
            model += ptr->eval_gradients(event, *this);

            // Optionally set gradient vector
            if (gradient != NULL) {
                for (int k = 0; k < ptr->size(); ++k, ++igrad) {
                    double grad        = (*ptr)(k).gradient();
                    (*gradient)(igrad) = (std::isinf(grad)) ? 0.0 : grad;
                }
            }
                
        } // endif: model was applicable

        // ... otherwise set gradient vector to 0
        else if (gradient != NULL) {
            for (int k = 0; k < ptr->size(); ++k, ++igrad) {
                (*gradient)(igrad) = 0.0;
            }
        }

    } // endfor: looped over models

    // Return
    return model;
}


/***********************************************************************//**
 * @brief Return total number of predicted counts for all models.
 *
 * @param[in] models Models.
 * @param[out] gradient Model parameter gradients.
 *
 * Returns the total number of predicted counts within the analysis region.
 ***************************************************************************/
double GObservation::npred(const GModels& models, GVector* gradient) const
{
    // Initialise
    double npred = 0.0;    // Reset predicted number of counts
    int    igrad = 0;      // Reset gradient counter

    // Loop over models
    for (int i = 0; i < models.size(); ++i) {

        // Extract pointer to model (circumvent const-correctness)
        GModel* model = (GModel*)models(i);

        // Handle only components that are relevant for the actual
        // instrument
        if (model->isvalid(instrument())) {

            // Determine Npred for model
            npred += npred_temp(*model);

            // Determine Npred gradients (perform computation only for free
            // parameters)
            for (int k = 0; k < model->size(); ++k) {
                if ((*model)(k).isfree()) {
                    double grad = npred_grad_temp(*model, k);
                    (*gradient)(igrad+k) = grad;
                }
            }

        } // endif: model component was valid for instrument

        // Increment parameter counter for gradient
        igrad += model->size();

    } // endfor: Looped over models

    // Return prediction
    return npred;
}


/***********************************************************************//**
 * @brief Set observation name
 *
 * @param[in] obsname Observation name.
 *
 * Set the name of this observation.
 ***************************************************************************/
void GObservation::obsname(const std::string& obsname)
{
    // Set name
    m_obsname = obsname;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set region of interest
 *
 * @param[in] roi Region of interest.
 *
 * Set region of interest for this observation.
 ***************************************************************************/
void GObservation::roi(const GRoi& roi)
{
    // Remove an existing ROI
    if (m_roi != NULL) delete m_roi;
    
    // Set ROI
    m_roi = roi;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set events
 *
 * @param[in] events Events.
 *
 * Set events for this observation.
 ***************************************************************************/
void GObservation::events(const GEvents& events)
{
    // Remove an existing events
    if (m_events != NULL) delete m_events;
    
    // Set events
    m_events = events;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set optimizer statistics
 *
 * @param[in] statistics Optimizer statistics.
 *
 * Set the optimizer statistics for this observation.
 ***************************************************************************/
void GObservation::statistics(const std::string& statistics)
{
    // Set statistics
    m_statistics = statistics;

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
    m_obsname.clear();
    m_ebounds.clear();
    m_gti.clear();
    m_roi        = NULL;
    m_events     = NULL;
    m_statistics = "Poisson";

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
    // Copy attributes
    m_obsname    = obs.m_obsname;
    m_ebounds    = obs.m_ebounds;
    m_gti        = obs.m_gti;
    m_statistics = obs.m_statistics;

    // Clone members that exist
    m_roi    = (obs.m_roi    != NULL) ? obs.m_roi->clone()      : NULL;
    m_events = (obs.m_events != NULL) ? obs.m_events->clone()   : NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GObservation::free_members(void)
{
    // Free memory
    if (m_events != NULL) delete m_events;
    if (m_roi    != NULL) delete m_roi;

    // Signal free pointers
    m_events = NULL;
    m_roi    = NULL;

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                        Npred integration methods                        =
 =                                                                         =
 ==========================================================================*/

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
 * \f[\int_{\rm GTI} \int_{E_{\rm bounds}} \int_{\rm ROI} S(\vec{p}, E, t)
 *    PSF(\vec{p'} | \vec{d}, \vec{p}, E, t) \, {\rm d}\vec{p'}\f]
 * where
 * \f$S(\vec{p}, E, t)\f$ is the source model,
 * \f$PSF(\vec{p'} | \vec{d}, \vec{p}, E, t) \, {\rm d}\vec{p'}\f$ is the
 * point spread function,
 * \f$\vec{p}\f$ is the true photon arrival direction,
 * \f$\vec{p'}\f$ is the measured photon direction,
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
    for (int i = 0; i < m_gti.size(); ++i) {

        // Set integration interval in MET
        double tstart = m_gti.tstart(i).met();
        double tstop  = m_gti.tstop(i).met();

        // Throw exception if time interval is not valid
        if (tstop <= tstart)
            throw GException::gti_invalid(G_NPRED_TEMP, m_gti.tstart(i),
                                          m_gti.tstop(i));

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
 * @param[in] srcTime True photon arrival time.
 *
 * @exception GException::erange_invalid
 *            Energy range is invalid.
 *
 * Computes
 * \f[\int_{E_{\rm bounds}} \int_{\rm ROI} S(\vec{p}, E, t)
 *    PSF(\vec{p'} | \vec{d}, \vec{p}, E, t) \, {\rm d}\vec{p'}\f]
 * where
 * \f$S(\vec{p}, E, t)\f$ is the source model,
 * \f$PSF(\vec{p'} | \vec{d}, \vec{p}, E, t) \, {\rm d}\vec{p'}\f$ is the
 * point spread function,
 * \f$\vec{p}\f$ is the true photon arrival direction,
 * \f$\vec{p'}\f$ is the measured photon direction,
 * \f$E\f$ is the true photon energy,
 * \f$t\f$ is the true photon arrival time, and
 * \f$d\f$ is the instrument pointing.
 *
 * \f$E_{\rm bounds}\f$ are the energy boundaries that are stored in the
 * GObservation::m_ebounds member.
 ***************************************************************************/
double GObservation::npred_spec(const GModel& model,
                                const GTime& srcTime) const
{
    // Set integration energy interval in MeV
    double emin = m_ebounds.emin().MeV();
    double emax = m_ebounds.emax().MeV();

    // Throw exception if energy range is not valid
    if (emax <= emin)
        throw GException::erange_invalid(G_NPRED_SPEC, emin, emax);

    // Setup integration function
    GObservation::npred_kern_spat integrand(this, model, srcTime);
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
 * @brief Spatially integrates the Npred kernel
 *
 * @param[in] model Gamma-ray source model.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 *
 * @exception GException::no_response
 *            No valid instrument response function for observation.
 * @exception GException::feature_not_implemented
 *            Computation for data-space models not yet implemented.
 *
 * Computes
 * \f[\int_{\rm ROI} S(\vec{p}, E, t)
 *    PSF(\vec{p'} | \vec{d}, \vec{p}, E, t) \, {\rm d}\vec{p'}\f]
 * where
 * \f$S(\vec{p}, E, t)\f$ is the source model,
 * \f$PSF(\vec{p'} | \vec{d}, \vec{p}, E, t) \, {\rm d}\vec{p'}\f$ is the
 * point spread function,
 * \f$\vec{p}\f$ is the true photon arrival direction,
 * \f$\vec{p'}\f$ is the measured photon direction,
 * \f$E\f$ is the true photon energy,
 * \f$t\f$ is the true photon arrival time, and
 * \f$d\f$ is the instrument pointing.
 *
 * \f${\rm ROI}\f$ is the region of interest that is stored in the
 * GObservation::m_roi member.
 * Note that the integration is performed by the GResponse::npred() method.
 *
 * @todo Implement computation for data-space model
 ***************************************************************************/
double GObservation::npred_spat(const GModel& model,
                                const GEnergy& srcEng,
                                const GTime& srcTime) const
{
    // Initialise result
    double result = 0.0;

    // Get sky model (NULL if not a sky model)
    GModelSky* sky = dynamic_cast<GModelSky*>((GModel*)&model);

    // Case A: We have a sky model
    if (sky != NULL) {

        // Get response function
        GResponse* rsp = response();
        if (rsp == NULL)
            throw GException::no_response(G_NPRED_SPAT);

        // Compute integrated IRF
        double npred = rsp->npred(*sky, srcEng, srcTime, *this);

        // Compute source model
        double source = 1.0;
        if (sky->spectral() != NULL) npred  *= sky->spectral()->eval(srcEng);
        if (sky->temporal() != NULL) source *= sky->temporal()->eval(srcTime);

        // Set result
        result = source * npred;
    
    }

    // Case B: We have a data space model
    else {
        throw GException::feature_not_implemented(G_NPRED_SPAT);
    }

    // Return
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
    double value = m_parent->npred_spat(*m_model, eng, *m_time);

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
 =                    Npred gradient integration methods                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Temporally integrate spatially/spectrally integrated Npred kernel
 *
 * @param[in] model Gamma-ray source model.
 * @param[in] ipar Parameter index for which gradient should be returned.
 *
 * @exception GException::gti_invalid
 *            Good Time Interval invalid.
 *
 * Note that MET is used for the time integration interval. This, however,
 * is no specialisation since npred_grad_kern_spec::eval() converts the
 * argument back in a GTime object by assuming that the argument is in MET,
 * hence the correct time system will be used at the end by the method.
 ***************************************************************************/
double GObservation::npred_grad_temp(const GModel& model, int ipar) const
{
    // Initialise result
    double result = 0.0;

    // Loop over GTIs
    for (int i = 0; i < m_gti.size(); ++i) {

        // Set integration interval in MET
        double tstart = m_gti.tstart(i).met();
        double tstop  = m_gti.tstop(i).met();

        // Throw exception if time interval is not valid
        if (m_gti.tstop(i) <= m_gti.tstart(i))
            throw GException::gti_invalid(G_NPRED_GRAD_TEMP, m_gti.tstart(i),
                                          m_gti.tstop(i));

        // Setup integration function
        GObservation::npred_grad_kern_spec integrand(this, model, ipar);
        GIntegral                          integral(&integrand);

        // Do Romberg integration
        result += integral.romb(tstart, tstop);

    } // endfor: looped over GTIs

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Spectrally integrate spatially integrated Npred gradient kernel
 *
 * @param[in] model Gamma-ray source model.
 * @param[in] ipar Parameter index for which gradient should be returned.
 * @param[in] srcTime True photon arrival time.
 *
 * @exception GException::erange_invalid
 *            Energy range is invalid.
 ***************************************************************************/
double GObservation::npred_grad_spec(const GModel& model, int ipar,
                                     const GTime& srcTime) const
{
    // Set integration energy interval in TeV
    double emin = m_ebounds.emin().MeV();
    double emax = m_ebounds.emax().MeV();

    // Throw exception if energy range is not valid
    if (emax <= emin)
        throw GException::erange_invalid(G_NPRED_GRAD_SPEC, emin, emax);

    // Setup integration function
    GObservation::npred_grad_kern_spat integrand(this, model, ipar, srcTime);
    GIntegral                          integral(&integrand);

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
 * @brief Integrates the Npred gradient kernel spatially
 *
 * @param[in] model Gamma-ray source model.
 * @param[in] ipar Parameter index for which gradient should be returned.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt Instrument pointing direction.
 *
 * @exception GException::no_response
 *            No valid instrument response function for observation.
 * @exception GException::feature_not_implemented
 *            Method not yet implemented for data-space model.
 *
 * @todo Implement method for data-space model.
 * @todo Avoid recomputation of gradients for same model parameters.
 ***************************************************************************/
double GObservation::npred_grad_spat(const GModel& model, int ipar,
                                     const GEnergy& srcEng,
                                     const GTime& srcTime) const
{
    // Initialise result
    double result = 0.0;

    // Get sky model (NULL if not a sky model)
    GModelSky* sky = dynamic_cast<GModelSky*>((GModel*)&model);

    // Case A: We have a sky model
    if (sky != NULL) {

        // Get response function
        GResponse* rsp = response();
        if (rsp == NULL)
            throw GException::no_response(G_NPRED_GRAD_SPAT);

        // Compute integrated IRF
        double npred = rsp->npred(*sky, srcEng, srcTime, *this);

        // Compute source model
        double source = 1.0;
        if (sky->spectral() != NULL) source *= sky->spectral()->eval_gradients(srcEng);
        if (sky->temporal() != NULL) source *= sky->temporal()->eval_gradients(srcTime);

        // Extract requested model gradient
        double grad = (*sky)(ipar).gradient();

        // Compute result
        result = grad * npred;

    }

    // Case B: We have a data space model
    else {
        throw GException::feature_not_implemented(G_NPRED_SPAT);
    }
    
    // Return
    return result;
}


/***********************************************************************//**
 * @brief Integration kernel for npred_grad_spec() method
 *
 * @param[in] x Function value.
 *
 * This method implements the integration kernel needed for the
 * npred_grad_spec() method. Upon the defintion of the G_LN_ENERGY_INT
 * declaration the energy integration is done logarithmically
 * (G_LN_ENERGY_INT=1) or not.
 ***************************************************************************/
double GObservation::npred_grad_kern_spat::eval(double x)
{
    #if G_LN_ENERGY_INT
    // Variable substitution
    x = exp(x);
    #endif

    // Set energy in MeV
    GEnergy eng;
    eng.MeV(x);

    // Get function value
    double value = m_parent->npred_grad_spat(*m_model, m_ipar, eng, *m_time);

    #if G_LN_ENERGY_INT
    // Correct for variable substitution
    value *= x;
    #endif

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Integration kernel for npred_grad_temp() method
 *
 * @param[in] x Function value.
 *
 * Note that MET is used for the time conversion. This, however, is no
 * specialisation since npred_grad_temp() hands MET.
 ***************************************************************************/
double GObservation::npred_grad_kern_spec::eval(double x)
{
    // Convert argument in MET
    GTime time;
    time.met(x);

    // Return value
    return (m_parent->npred_grad_spec(*m_model, m_ipar, time));
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
