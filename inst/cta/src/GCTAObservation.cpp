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
#include "GException.hpp"
#include "GCTAObservation.hpp"
#include "GCTAEventList.hpp"
#include "GCTAEventCube.hpp"
#include "GCTARoi.hpp"
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
 * @param[in] filename Event FITS file name.
 ***************************************************************************/
void GCTAObservation::load_unbinned(const std::string& filename)
{
    // Delete old events
    if (m_events != NULL) delete m_events;

    // Allocate event list and establish link to observation
    m_events = new GCTAEventList;
    m_events->obs(this);

    // Load events into list
    m_events->load(filename);

    // Load GTIs
    //m_gti.load(filename);

    // Set attributes
    //m_tstart.met(m_gti.tstart());
    //m_tstop.met(m_gti.tstop());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load data for binned analysis
 *
 * @param[in] filename Counts map FITS file name.
 ***************************************************************************/
void GCTAObservation::load_binned(const std::string& filename)
{
    // Delete old events
    if (m_events != NULL) delete m_events;

    // Allocate event cube and establish link to observation
    m_events = new GCTAEventCube;
    m_events->obs(this);

    // Load events into cube
    m_events->load(filename);

    // Copy energy boundaries and GTIs from event cube
    m_ebounds = ((GCTAEventCube*)m_events)->m_ebds;
    m_gti     = ((GCTAEventCube*)m_events)->m_gti;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return total number of predicted counts for all models.
 *
 * @param[in] models Models.
 * @param[in] gradient Model parameter gradients.
 ***************************************************************************/
double GCTAObservation::npred(const GModels& models, GVector* gradient) const
{
    // Initialise
    double npred = 0.0;    // Reset predicted number of counts
    int    igrad = 0;      // Reset gradient counter

    // Loop over models
    for (int i = 0; i < models.size(); ++i) {

        // Extract pointer to model (bypass const-correctness)
        GModel* model = (GModel*)models(i);

        // Handle only components that are relevant for CTA
        if (model->isvalid("CTA")) {

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


/*==========================================================================
 =                                                                         =
 =                        Npred integration methods                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Temporally integrate spatially & spectrally integrated Npred kernel
 *
 * @param[in] model Gamma-ray source model.
 *
 * Implement the temporal integration as a simple multiplication by the
 * elapsed time. This assumes that the source is non-variable during the
 * observation and that the CTA pointing is stable.
 ***************************************************************************/
double GCTAObservation::npred_temp(const GModel& model) const
{
    // Initialise result
    double result = 0.0;

    // Determine ontime
    double ontime = m_gti.ontime();
    
    // Integrate only if ontime is positive
    if (ontime > 0.0) {

        // Integration is a simple multiplication by the time
        result = npred_spec(model, m_gti.tstart()) * ontime;

    }

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                    Npred gradient integration methods                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Temporally integrate spatially & spectrally integrated Npred gradient kernel
 *
 * @param[in] model Gamma-ray source model.
 * @param[in] ipar Parameter index for which gradient should be returned.
 *
 * Implement the temporal integration as a simple multiplication by the
 * elapsed time. This assumes that the source is non-variable during the
 * observation and that the CTA pointing is stable.
 ***************************************************************************/
double GCTAObservation::npred_grad_temp(const GModel& model, int ipar) const
{
    // Initialise result
    double result = 0.0;

    // Determine ontime
    double ontime = m_gti.ontime();
    
    // Integrate only if ontime is positive
    if (ontime > 0.0) {

        // Integration is a simple multiplication by the time
        result = npred_grad_spec(model, ipar, m_gti.tstart()) * ontime;

    }

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
    os << " Time range ................: " << std::fixed
       << obs.m_gti.tstart().mjd() << " - "
       << obs.m_gti.tstop().mjd() << " days" << std::endl;
    os << " Energy range ..............: " << std::fixed
       << obs.m_ebounds.emin().MeV() << " - "
       << obs.m_ebounds.emax().MeV() << " MeV" << std::endl;

    // Add ROI to stream if it exists
    if (obs.m_roi != NULL)
        os << " Region of interest ........: " << *((GCTARoi*)obs.m_roi) << std::endl;

    // Add GTIs to stream
    os << obs.m_gti << std::endl;

    // Add response to stream if it exists
    if (obs.m_response != NULL)
        os << *((GCTAResponse*)obs.m_response) << std::endl;

    // Add events to stream
    if (obs.m_events != NULL) {
        if (obs.m_events->islist())
            os << *((GCTAEventList*)obs.m_events) << std::endl;
        else
            os << *((GCTAEventCube*)obs.m_events) << std::endl;
    }

    // Return output stream
    return os;
}
