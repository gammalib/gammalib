/***************************************************************************
 *               GCTAObservation.cpp  -  CTA Observation class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GCTAObservation.cpp
 * @brief GCTAObservation class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include "GCTAException.hpp"
#include "GCTAObservation.hpp"
#include "GCTAEventList.hpp"
#include "GFits.hpp"
#include "GModelSpatialPtsrc.hpp"
#include "GIntegral.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 *
 * Creates an empty instance of GCTAObservation.
 ***************************************************************************/
GCTAObservation::GCTAObservation(void) : GObservation()
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
 *
 * Creates an instance of GCTAObservation by copying information from an
 * existing instance.
 ***************************************************************************/
GCTAObservation::GCTAObservation(const GCTAObservation& obs) : GObservation(obs)
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
GCTAObservation::~GCTAObservation(void)
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
 *
 * Assign instance of GCTAObservation to this object.
 ***************************************************************************/
GCTAObservation& GCTAObservation::operator= (const GCTAObservation& obs)
{
    // Execute only if object is not identical
    if (this != &obs) {

        // Copy base class members
        this->GObservation::operator=(obs);

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
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Set CTA response function
 *
 * @param[in] irfname Name of CTA response function.
 * @param[in] caldb Optional name of calibration database.
 *
 * @todo Add a generic method in GResponse that searches for path to
 *       calibration database using environment variables.
 ***************************************************************************/
void GCTAObservation::response(const std::string& irfname, std::string caldb)
{
    // Delete old response function
    if (m_response != NULL) delete m_response;
    
    // Allocate new CTA response function
    m_response = new GCTAResponse;
    
    // Set calibration database
    m_response->caldb(caldb);
    
    // Load instrument response function
    m_response->load(irfname);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Load data for unbinned analysis
 *
 * @param[in] evname Event FITS filename.
 ***************************************************************************/
void GCTAObservation::load_unbinned(const std::string& evname)
{
    // Delete old events
    if (m_events != NULL) delete m_events;

    // Allocate events and establish link to observation
    m_events = new GCTAEventList;
    m_events->obs(this);
    
    // Load events
    m_events->load(evname);

    // Load GTIs
    //m_gti.load(evname);

    // Set attributes
    //m_tstart.met(m_gti.tstart());
    //m_tstop.met(m_gti.tstop());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return total number of predicted counts for all models.
 *
 * @param[in] models Models.
 *
 * @todo Method needs to be implemented
 ***************************************************************************/
double GCTAObservation::npred(const GModels& models) const
{
    // Initialise predicted number of counts
    double npred = 0.0;
    
    // Loop over models
    for (int i = 0; i < models.size(); ++i) {
    
        // Check if model is a CTA model and if it should be used for this
        // event
        // TO BE IMPLEMENTED
        
        // Add model
        GCTAPointing pnt;
        std::cout << *models(i) << std::endl;
        npred += npred_integrate_temporal(*models(i), pnt);
        
    }

    // Return prediction
    return npred;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCTAObservation::init_members(void)
{
    // Set instrument name
    m_instrument = "CTA";

    // Initialise members

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] obs Observation to be copied
 ***************************************************************************/
void GCTAObservation::copy_members(const GCTAObservation& obs)
{
    // Copy members

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAObservation::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone class
***************************************************************************/
GCTAObservation* GCTAObservation::clone(void) const
{
    return new GCTAObservation(*this);
}


/***********************************************************************//**
 * @brief Computes the Npred integrand
 *
 * @param[in] model Gamma-ray source model.
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt Instrument pointing direction.
 *
 * Computes the integrand for the Npred computation. This method is called
 * by npred_integrate_spatial() which performs the spatial integration of
 * integrand.
 ***************************************************************************/
double GCTAObservation::npred_integrand(const GModel& model,
                                        const GSkyDir& srcDir,
                                        const GEnergy& srcEng,
                                        const GTime& srcTime,
                                        const GPointing& pnt) const
{
    // Compute integrated IRF
    double nirf = m_response->nirf(srcDir, srcEng, srcTime, pnt);

    // Compute source model
    GCTAInstDir obsDir;  // unused, to be removed later
    GEnergy     obsEng;  // unused, to be removed later
    GTime       obsTime; // unused, to be removed later
    GModel* ptr    = (GModel*)&model; // impose const-correctness
    double  source = ptr->value(obsDir, obsEng, obsTime, srcDir, srcEng,
                                srcTime, *m_response, pnt);

    // Return
    return (source * nirf);
}


/***********************************************************************//**
 * @brief Integrates the Npred integrand spatially
 *
 * @param[in] model Gamma-ray source model.
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt Instrument pointing direction.
 ***************************************************************************/
double GCTAObservation::npred_integrate_spatial(const GModel& model,
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
        std::cout << "GCTAObservation::npred_integrate_spatial:"
                  << " Integration not implemented." << std::endl;
    }
    
    // Case B: No integration, then extract point source position from model
    else {        
        GSkyDir srcDir;
        srcDir.radec_deg(((GModelSpatialPtsrc*)model.spatial())->ra(),
                         ((GModelSpatialPtsrc*)model.spatial())->dec());
        result = npred_integrand(model, srcDir, srcEng, srcTime, pnt);
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Integrates the Npred integrand spectrally
 *
 * @param[in] model Gamma-ray source model.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt Instrument pointing direction.
 ***************************************************************************/
double GCTAObservation::npred_integrate_spectral(const GModel& model,
                                                 const GTime& srcTime,
                                                 const GPointing& pnt) const
{
    // Setup integration function
    GCTAObservation::int_spec integrand(this, model, srcTime, pnt);
    GIntegral                 integral(&integrand);
    
    // Do Romberg integration
    double result = integral.romb(0.02, 100.0); // energies in TeV
    
    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Integrates the Npred integrand temporally
 *
 * @param[in] model Gamma-ray source model.
 * @param[in] pnt Instrument pointing direction.
 ***************************************************************************/
double GCTAObservation::npred_integrate_temporal(const GModel& model,
                                                 const GPointing& pnt) const
{
    // Setup integration function
    GCTAObservation::int_temp integrand(this, model, pnt);
    GIntegral                 integral(&integrand);
    
    // Do Romberg integration
    double result = integral.romb(0.0, 100000.0); // 100000 sec
    
    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put CTA observation in output stream
 *
 * @param[in] os Output stream into which the data will be dumped
 * @param[in] obs Observation to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GCTAObservation& obs)
{
    // Put observation in stream
    os.precision(3);
    os << "=== GCTAObservation ===" << std::endl;
    os << " Name ......................: " << obs.m_obsname << std::endl;
    os << " Instrument ................: " << obs.m_instrument << std::endl;
    os << " Time range ................: " << std::fixed << 
            obs.m_tstart.mjd() << " - " << obs.m_tstop.mjd() << std::endl;
    os << " Energy range ..............: " << std::fixed << 
            obs.m_emin.TeV() << " - " << obs.m_emax.TeV() << " TeV" << std::endl;

    // Add events to stream
    if (obs.m_events != NULL) {
        if (obs.m_events->islist())
            os << *((GCTAEventList*)obs.m_events) << std::endl;
        else
            os << std::endl;
//            os << *((GCTAEventCube*)obs.m_events) << std::endl;
    }

    // Add GTIs to stream
    os << obs.m_gti << std::endl;
    
    // Return output stream
    return os;
}
