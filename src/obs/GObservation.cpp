/***************************************************************************
 *           GObservation.cpp  -  Observation abstract base class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2010 by Jurgen Knodlseder                           *
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
#include <iostream>
#include <cmath>
#include "GException.hpp"
#include "GObservation.hpp"
#include "GModelSpatialPtsrc.hpp"
#include "GIntegral.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_MODEL           "GObservation::model(GModels&, GPointing&, GInstDir&, GEnergy&, GTime&, GVector*) const"
#define G_NPRED_KERN      "GObservation::npred_kern(GModel&, GSkyDir&, GEnergy&, GTime&, GPointing&)"
#define G_NPRED_SPEC      "GObservation::npred_spec(GModel&, GTime&)"
#define G_NPRED_TEMP      "GObservation::npred_temp(GModel&, int)"
#define G_NPRED_GRAD_KERN "GObservation::npred_grad_kern(GModel&, int, GSkyDir&, GEnergy&, GTime&, GPointing&)"
#define G_NPRED_GRAD_SPEC "GObservation::npred_grad_spec(GModel&, int, GTime&)"
#define G_NPRED_GRAD_TEMP "GObservation::npred_grad_temp(GModel&, int)"

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
 * @param[in] obsDir Measured instrument direction.
 * @param[in] obsEng Measured energy.
 * @param[in] obsTime Measured arrival time.
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
double GObservation::model(const GModels& models, const GInstDir& obsDir,
                           const GEnergy& obsEng, const GTime& obsTime,
                           GVector* gradient) const
{
    // Verify that gradients vector has the same dimension than the
    // model has parameters
    #if defined(G_RANGE_CHECK)
    if (models.npars() != gradient->size())
        throw GException::gradient_par_mismatch(G_MODEL, gradient->size(),
                                                models.npars());
    #endif

    // Initialise model value
    double model = 0.0;

    // Get current pointing and response function
    GPointing* pnt = pointing(obsTime);
    GResponse* rsp = response(obsTime);

    // Loop over models
    for (int i = 0; i < models.size(); ++i) {

        // Get model pointer (enforce const correctness)
        GModel* ptr = (GModel*)models(i);

        // Check if model applies to specific instrument
        if (ptr->isvalid(instrument())) {

            // Compute value and add to model
            model += ptr->eval_gradients(obsDir, obsEng, obsTime, *rsp, *pnt);

        } // endif: model was applicable

    } // endfor: looped over models

    // Optionally set gradient vector
    if (gradient != NULL) {
        for (int i = 0; i < gradient->size(); ++i) {
            double grad    = models.par(i)->gradient();
            (*gradient)(i) = (std::isinf(grad)) ? 0.0 : grad;
        }
    }

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

        // Extract pointer to model (bypass const-correctness)
        GModel* model = (GModel*)models(i);

        // Handle only components that are relevant for the actual
        // instrument
        if (model->isvalid(instrument())) {

            // Determine Npred for model
            npred += npred_temp(*model);

            // Determine Npred gradients (perform computation only for free
            // parameters)
            for (int k = 0; k < model->npars(); ++k) {
                if (model->par(k)->isfree()) {
                    double grad = npred_grad_temp(*model, k);
                    (*gradient)(igrad+k) = grad;
                }
            }

        } // endif: model component was valid for instrument

        // Increment parameter counter for gradient
        igrad += model->npars();

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
    m_roi    = (obs.m_roi      != NULL) ? obs.m_roi->clone()      : NULL;
    m_events = (obs.m_events   != NULL) ? obs.m_events->clone()   : NULL;

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
 * @brief Temporally integrates spatially and spectrally integrated Npred kernel
 *
 * @param[in] model Gamma-ray source model.
 *
 * @exception GException::gti_invalid
 *            Good Time Interval is invalid.
 *
 *
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
 * Integrates the Npred kernel over the sky region of interest. In case that
 * the sky model is a point-source no integration is performed and the
 * kernel value is directly returned.
 *
 * @todo Implement integration over skymap.
 ***************************************************************************/
double GObservation::npred_spat(const GModel& model,
                                const GEnergy& srcEng,
                                const GTime& srcTime) const
{
    // Initialise result
    double result = 0.0;

    // Continue only if the gamma-ray source model has a spatial component
    if (model.spatial() != NULL) {

        // Case A: Model is a point source
        if (model.spatial()->isptsource()) {

            // Build sky direction from point source parameters
            GSkyDir srcDir;
            srcDir.radec_deg(((GModelSpatialPtsrc*)model.spatial())->ra(),
                             ((GModelSpatialPtsrc*)model.spatial())->dec());

            // Get function value at that position
            result = npred_kern(model, srcDir, srcEng, srcTime);

        } // endif: Model was a point source

        // Case B: Model is not a point source
        else {

            // Dump warning that integration is not yet implemented
            std::cout << "WARNING: GObservation::npred_spat:"
                      << " Sky integration not implemented." << std::endl;

        } // endelse: Model was not a point source

    } // endif: Gamma-ray source model had a spatial component

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Computes the Npred kernel
 *
 * @param[in] model Gamma-ray source model.
 * @param[in] srcDir True photon arrival direction.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 *
 * @exception GException::roi_invalid
 *            No valid ROI found for observation.
 *
 * Computes the integrand (or integral kernel) for the Npred computation.
 * This method is called by npred_spat() which performs the spatial
 * integration of the kernel.
 ***************************************************************************/
double GObservation::npred_kern(const GModel& model,
                                const GSkyDir& srcDir,
                                const GEnergy& srcEng,
                                const GTime& srcTime) const
{
    // Make sure that ROI has been set
    if (m_roi == NULL)
        throw GException::roi_invalid(G_NPRED_KERN,
                          "No ROI has been defined for observation.");

    // Get current pointing and response
    GPointing* pnt = pointing(srcTime);
    GResponse* rsp = response(srcTime);

    // Compute integrated IRF
    double nirf = rsp->nirf(srcDir, srcEng, srcTime, *pnt, *m_roi, m_ebounds, m_gti);

    // Compute source model
    GModel* ptr    = (GModel*)&model; // bypass const-correctness
    double  source = ptr->value(srcDir, srcEng, srcTime);

    // Return
    return (source * nirf);
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
 * @todo Throw exception if integration energy range is not valid.
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
 * @todo Implement integration over skymap.
 ***************************************************************************/
double GObservation::npred_grad_spat(const GModel& model, int ipar,
                                     const GEnergy& srcEng,
                                     const GTime& srcTime) const
{
    // Initialise result
    double result = 0.0;

    // Continue only if the gamma-ray source model has a spatial component
    if (model.spatial() != NULL) {

        // Case A: Model is a point source
        if (model.spatial()->isptsource()) {

            // Build sky direction from point source parameters
            GSkyDir srcDir;
            srcDir.radec_deg(((GModelSpatialPtsrc*)model.spatial())->ra(),
                             ((GModelSpatialPtsrc*)model.spatial())->dec());

            // Get function value at that position
            result = npred_grad_kern(model, ipar, srcDir, srcEng, srcTime);

        } // endif: Model was a point source

        // Case B: Model is not a point source
        else {

            // Dump warning that integration is not yet implemented
            std::cout << "WARNING: GObservation::npred_grad_spat:"
                      << " Sky integration not implemented." << std::endl;

        } // endelse: Model was not a point source

    } // endif: Gamma-ray source model had a spatial component

    // Return result
    return result;
}



/***********************************************************************//**
 * @brief Returns the Npred gradient for a given model parameter
 *
 * @param[in] model Gamma-ray source model.
 * @param[in] ipar Parameter index for which gradient should be returned.
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 *
 * @exception GException::roi_invalid
 *            No valid ROI found for observation.
 ***************************************************************************/
double GObservation::npred_grad_kern(const GModel& model, int ipar,
                                     const GSkyDir& srcDir,
                                     const GEnergy& srcEng,
                                     const GTime& srcTime) const
{
    // Make sure that ROI has been set
    if (m_roi == NULL)
        throw GException::roi_invalid(G_NPRED_GRAD_KERN,
                          "No ROI has been defined for observation.");

    // Get current pointing and response
    GPointing* pnt = pointing(srcTime);
    GResponse* rsp = response(srcTime);

    // Compute integrated IRF
    double nirf = rsp->nirf(srcDir, srcEng, srcTime, *pnt, *m_roi, m_ebounds, m_gti);

    // Get model gradients
    GModel* ptr       = (GModel*)&model; // bypass const-correctness
    GVector gradients = ptr->gradients(srcDir, srcEng, srcTime);

    // Extract requested model gradient
    double grad = gradients(ipar);

    // Return
    return (grad * nirf);
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
 * @param[in] os Output stream into which the data will be dumped
 * @param[in] obs Observation to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GObservation& obs)
{
    // Put observation in stream
    os.precision(3);
    os << "=== GObservation ===" << std::endl;
    os << " Name ......................: " << obs.m_obsname << std::endl;
    os << " Instrument ................: " << obs.instrument() << std::endl;
    os << " Time range ................: " << std::fixed
       << obs.m_gti.tstart().mjd() << " - "
       << obs.m_gti.tstop().mjd() << " days" << std::endl;
    os << " Energy range ..............: " << std::fixed
       << obs.m_ebounds.emin().MeV() << " - "
       << obs.m_ebounds.emax().MeV() << " MeV" << std::endl;
    os << " Optimizer statistics ......: " << obs.m_statistics << std::endl;

    // Add event list to stream
    if (obs.m_events != NULL)
        os << *(obs.m_events) << std::endl;

    // Add energy intervals to stream
    os << obs.m_ebounds << std::endl;

    // Add GTIs to stream
    os << obs.m_gti << std::endl;

    // Return output stream
    return os;
}
