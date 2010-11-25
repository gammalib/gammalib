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
#include "GCTAEventCube.hpp"
#include "GCTARoi.hpp"
#include "GException.hpp"
#include "GFits.hpp"
#include "GModelSpatialPtsrc.hpp"
#include "GIntegral.hpp"
#include "GIntegrand.hpp"

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
 * @brief Clone object
***************************************************************************/
GCTAObservation* GCTAObservation::clone(void) const
{
    return new GCTAObservation(*this);
}


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
 * @brief Returns pointer to CTA response function
 *
 * @param[in] time Time.
 *
 * Returns pointer to response function for a given time. As the response is
 * supposed not to vary during an observation, the time argument needs not to
 * be considered.
 ***************************************************************************/
GResponse* GCTAObservation::response(const GTime& time) const
{
    // Return response pointer
    return m_response;
}


/***********************************************************************//**
 * @brief Returns pointer to CTA pointing direction
 *
 * @param[in] time Time.
 *
 * Returns pointer to pointing direction for a given time. As the pointing
 * direction is supposed not to vary during an observation, the time argument
 * needs not to be considered.
 ***************************************************************************/
GPointing* GCTAObservation::pointing(const GTime& time) const
{
    // Return response pointer
    return m_pointing;
}


/***********************************************************************//**
 * @brief Returns instrument name
 ***************************************************************************/
std::string GCTAObservation::instrument(void) const
{
    // Return instument name
    return ("CTA");
}


/***********************************************************************//**
 * @brief Load data for unbinned analysis
 *
 * @param[in] filename Event FITS file name.
 *
 * @todo Implement GTI loading.
 ***************************************************************************/
void GCTAObservation::load_unbinned(const std::string& filename)
{
    // Delete old events
    if (m_events != NULL) delete m_events;

    // Allocate events
    GCTAEventList* events = new GCTAEventList;
    m_events = events;

    // Load events into list
    events->load(filename);

    // Load GTIs
    //m_gti.load(filename);

    // Set attributes
    //m_tstart.met(m_gti.tstart());
    //m_tstop.met(m_gti.tstop());

    // Link observations to events. This has to be done after loading since
    // loading initialises the GCTAEventList object, hence resets the pointer
    // to the observation.
    //events->obs(this);

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

    // Allocate events
    GCTAEventCube* events = new GCTAEventCube;
    m_events = events;

    // Load events into cube
    events->load(filename);

    // Copy energy boundaries and GTIs from event cube
    m_ebounds = ((GCTAEventCube*)m_events)->m_ebds;
    m_gti     = ((GCTAEventCube*)m_events)->m_gti;

    // Link observations to events. This has to be done after loading since
    // loading initialises the GCTAEventCube object, hence resets the pointer
    // to the observation.
    //events->obs(this);

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
void GCTAObservation::init_members(void)
{
    // Initialise members
    m_response = NULL;
    m_pointing = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] obs Observation to be copied
 *
 * @todo Try to avoid the back pointer if possible, or flesh out a way that
 * makes this more solid. The point is: when events are copied the back
 * pointer to the observation needs to be updated, yet when the copying is
 * done in the events class the class does not known about the observation.
 * Thus, back pointer update has to be done by the observation class.
 ***************************************************************************/
void GCTAObservation::copy_members(const GCTAObservation& obs)
{
    // Copy members
    if (obs.m_response != NULL) m_response = obs.m_response->clone();
    if (obs.m_pointing != NULL) m_pointing = obs.m_pointing->clone();

    // Update the back pointer to link observation the actual observation
    // to the events. This has to be done here since the events that were
    // copied do not yet know to which observation they belong.
    /*
    if (m_events->islist())
        ((GCTAEventList*)m_events)->obs(this);
    else
        ((GCTAEventCube*)m_events)->obs(this);
    */

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAObservation::free_members(void)
{
    // Free memory
    if (m_response != NULL) delete m_response;
    if (m_pointing != NULL) delete m_pointing;

    // Mark memory as free
    m_response = NULL;
    m_pointing = NULL;

    // Return
    return;
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
    os << " Instrument ................: " << obs.instrument() << std::endl;
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
        os << *(obs.m_response) << std::endl;

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
