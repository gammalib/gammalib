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
#include "GException.hpp"
#include "GObservation.hpp"
#include "GModelSpatialPtsrc.hpp"
#include "GIntegral.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_NPRED_TEMP                  "GObservation::npred_temp(GModel&,int)"
#define G_NPRED_GRAD_TEMP        "GObservation::npred_grad_temp(GModel&,int)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
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
 * @param[in] obs Observation from which the instance should be built.
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
 * @param[in] obs Observation which should be assigned.
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
    m_instrument.clear();
    m_tstart.clear();
    m_tstop.clear();
    m_emin.clear();
    m_emax.clear();
    m_events   = NULL;
    m_response = NULL;
    m_gti      = GGti();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] obs GObservation members which should be copied.
 ***************************************************************************/
void GObservation::copy_members(const GObservation& obs)
{
    // Copy attributes
    m_obsname    = obs.m_obsname;
    m_instrument = obs.m_instrument;
    m_tstart     = obs.m_tstart;
    m_tstop      = obs.m_tstop;
    m_emin       = obs.m_emin;
    m_emax       = obs.m_emax;
    m_gti        = obs.m_gti;

    // Clone members that exist
    m_events   = (obs.m_events   != NULL) ? obs.m_events->clone()   : NULL;
    m_response = (obs.m_response != NULL) ? obs.m_response->clone() : NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GObservation::free_members(void)
{
    // Free memory
    if (m_events   != NULL) delete m_events;
    if (m_response != NULL) delete m_response;

    // Signal free pointers
    m_events   = NULL;
    m_response = NULL;

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                        Npred integration methods                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Computes the Npred kernel for a specific direction, energy and time
 *
 * @param[in] model Gamma-ray source model.
 * @param[in] srcDir Photon arrival direction.
 * @param[in] srcEng Photon energy.
 * @param[in] srcTime Photon arrival time.
 * @param[in] pnt Instrument pointing direction.
 *
 * Computes the integrand (or integral kernel) for the Npred computation.
 * This method is called by npred_spat() which performs the spatial
 * integration of the kernel.
 ***************************************************************************/
double GObservation::npred_kern(const GModel& model, const GSkyDir& srcDir,
                                const GEnergy& srcEng, const GTime& srcTime,
                                const GPointing& pnt) const
{
    // Compute integrated IRF
    double nirf = m_response->nirf(srcDir, srcEng, srcTime, pnt);

    // Compute source model
    GModel* ptr    = (GModel*)&model; // bypass const-correctness
    double  source = ptr->value(srcDir, srcEng, srcTime);

    // Return
    return (source * nirf);
}


/***********************************************************************//**
 * @brief Integrates the Npred kernel spatially
 *
 * @param[in] model Gamma-ray source model.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt Instrument pointing direction.
 *
 * @todo More save handling of point source model.
 ***************************************************************************/
double GObservation::npred_spat(const GModel& model, const GEnergy& srcEng,
                                const GTime& srcTime,
                                const GPointing& pnt) const
{
    // Initialise result
    double result = 0.0;

    // Determine if integration is needed
    bool integrate  = (model.spatial() != NULL) ? model.spatial()->depdir() : false;

    // Case A: Integraion
    if (integrate) {
        std::cout << "GObservation::npred_spat:"
                  << " Integration not implemented." << std::endl;
    }

    // Case B: No integration, then extract point source position from model
    else {
        GSkyDir srcDir;
        srcDir.radec_deg(((GModelSpatialPtsrc*)model.spatial())->ra(),
                         ((GModelSpatialPtsrc*)model.spatial())->dec());
        result = npred_kern(model, srcDir, srcEng, srcTime, pnt);
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Integrates spatially integrated  Npred kernel spectrally
 *
 * @param[in] model Gamma-ray source model.
 * @param[in] srcTime True photon arrival time.
 *
 * @todo Throw exception if integration energy range is not valid.
 * @todo Implement correct time dependent extraction of telescope
 *       pointing.
 ***************************************************************************/
double GObservation::npred_spec(const GModel& model, const GTime& srcTime) const
{
    // Get telescope pointing for a given time
    GPointing* pnt;

    // Setup integration function
    GObservation::npred_kern_spat integrand(this, model, srcTime, pnt);
    GIntegral                     integral(&integrand);

    // Set integration energy interval in TeV
    double emin = m_emin.TeV();
    double emax = m_emax.TeV();

    // Test validity of integration energy range
    //TODO

    // Do Romberg integration
    double result = integral.romb(emin, emax);

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Temporally integrates spatially and spectrally integrated Npred kernel
 *
 * @param[in] model Gamma-ray source model.
 *
 * @exception GException::gti_invalid
 *            Good Time Interval is invalid.
 ***************************************************************************/
double GObservation::npred_temp(const GModel& model) const
{
    // Set integration interval in MET
    double tstart = m_tstart.met();
    double tstop  = m_tstop.met();

    // Throw exception if time interval is not valid
    if (tstop <= tstart)
        throw GException::gti_invalid(G_NPRED_TEMP, &m_gti);


    // Setup integration function
    GObservation::npred_kern_spec integrand(this, model);
    GIntegral                     integral(&integrand);

    // Do Romberg integration
    double result = integral.romb(tstart, tstop);

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                    Npred gradient integration methods                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Returns the Npred gradient for a given model parameter
 *
 * @param[in] model Gamma-ray source model.
 * @param[in] ipar Parameter index for which gradient should be returned.
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt Instrument pointing direction.
 ***************************************************************************/
double GObservation::npred_grad_kern(const GModel& model, int ipar,
                                     const GSkyDir& srcDir,
                                     const GEnergy& srcEng,
                                     const GTime& srcTime,
                                     const GPointing& pnt) const
{
    // Compute integrated IRF
    double nirf = m_response->nirf(srcDir, srcEng, srcTime, pnt);

    // Get model gradients
    GModel* ptr       = (GModel*)&model; // bypass const-correctness
    GVector gradients = ptr->gradients(srcDir, srcEng, srcTime);

    // Extract requested model gradient
    double grad = gradients(ipar);

    // Return
    return (grad * nirf);
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
 * @todo More save handling of point source model.
 ***************************************************************************/
double GObservation::npred_grad_spat(const GModel& model, int ipar,
                                     const GEnergy& srcEng,
                                     const GTime& srcTime,
                                     const GPointing& pnt) const
{
    // Initialise result
    double result = 0.0;

    // Determine if integration is needed
    bool integrate  = (model.spatial() != NULL) ? model.spatial()->depdir() : false;

    // Case A: Integraion
    if (integrate) {
        std::cout << "GObservation::npred_grad_spat:"
                  << " Integration not implemented." << std::endl;
    }

    // Case B: No integration, then extract point source position from model
    else {
        GSkyDir srcDir;
        srcDir.radec_deg(((GModelSpatialPtsrc*)model.spatial())->ra(),
                         ((GModelSpatialPtsrc*)model.spatial())->dec());
        result = npred_grad_kern(model, ipar, srcDir, srcEng, srcTime, pnt);
    }

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
 * @todo Implement correct time dependent extraction of telescope
 *       pointing.
 ***************************************************************************/
double GObservation::npred_grad_spec(const GModel& model, int ipar,
                                     const GTime& srcTime) const
{
    // Set integration energy interval in TeV
    double emin = m_emin.TeV();
    double emax = m_emax.TeV();

    // Test validity of integration energy range
    //TODO

    // Set telescope pointing for given time
    GPointing *pnt;

    // Setup integration function
    GObservation::npred_grad_kern_spat integrand(this, model, ipar, srcTime, pnt);
    GIntegral                          integral(&integrand);

    // Do Romberg integration
    double result = integral.romb(emin, emax);

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Temporally integrate spatially/spectrally integrated Npred gradient kernel
 *
 * @param[in] model Gamma-ray source model.
 * @param[in] ipar Parameter index for which gradient should be returned.
 *
 * @exception GException::gti_invalid
 *            Good Time Interval invalid.
 ***************************************************************************/
double GObservation::npred_grad_temp(const GModel& model, int ipar) const
{
    // Set integration interval in MET
    double tstart = m_tstart.met();
    double tstop  = m_tstop.met();

    // Throw exception if time interval is not valid
    if (tstop <= tstart)
        throw GException::gti_invalid(G_NPRED_GRAD_TEMP, &m_gti);

    // Setup integration function
    GObservation::npred_grad_kern_spec integrand(this, model, ipar);
    GIntegral                          integral(&integrand);

    // Do Romberg integration
    double result = integral.romb(tstart, tstop);

    // Return result
    return result;
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
    os << " Instrument ................: " << obs.m_instrument << std::endl;
    os << " Time range ................: " << std::fixed << 
            obs.m_tstart.mjd() << " - " << obs.m_tstop.mjd() << std::endl;
    os << " Energy range ..............: " << std::fixed << 
            obs.m_emin.MeV() << " - " << obs.m_emax.MeV() << " MeV" << std::endl;

    // Add event list to stream
    if (obs.m_events != NULL)
        os << *(obs.m_events) << std::endl;

    // Add GTIs to stream
    os << obs.m_gti << std::endl;
    
    // Return output stream
    return os;
}
