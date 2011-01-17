/***************************************************************************
 *               GCTAObservation.cpp  -  CTA Observation class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
#include "GCTAException.hpp"
#include "GCTAObservation.hpp"
#include "GCTAEventList.hpp"
#include "GCTAEventCube.hpp"
#include "GCTARoi.hpp"
#include "GException.hpp"
#include "GFits.hpp"
#include "GTools.hpp"
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
 * @brief Clear instance
 ***************************************************************************/
void GCTAObservation::clear(void)
{
    // Free members
    free_members();
    this->GObservation::free_members();

    // Initialise members
    this->GObservation::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
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
 ***************************************************************************/
GCTAResponse* GCTAObservation::response(void) const
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
GCTAPointing* GCTAObservation::pointing(const GTime& time) const
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
 * @brief Print CTA observation information
 ***************************************************************************/
std::string GCTAObservation::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GCTAObservation ===\n");
    result.append(parformat("Name")+obsname()+"\n");
    result.append(parformat("Instrument")+instrument()+"\n");
    result.append(parformat("Statistics")+statistics()+"\n");

    // Append time range
    result.append(parformat("Time range"));
    result.append(str(m_gti.tstart().mjd()));
    result.append(" - ");
    result.append(str(m_gti.tstop().mjd()));
    result.append(" days\n");

    // Append energy range
    result.append(parformat("Energy range"));
    result.append(m_ebounds.emin().print());
    result.append(" - ");
    result.append(m_ebounds.emax().print());

    // Append ROI
    if (m_roi != NULL)
        result.append("\n"+m_roi->print());
    else
        result.append("\n"+parformat("Region of interest")+"undefined");

    // Append GTIs
    //result.append("\n"+m_gti->print());

    // Append response
    if (m_response != NULL)
        result.append("\n"+m_response->print());
    else
        result.append("\n"+parformat("CTA response")+"undefined");

    // Append events
    if (m_events != NULL)
        result.append("\n"+m_events->print());

    // Return result
    return result;
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
    // Delete old events. We do not call clear() here since we want to
    // preserve any existing response function.
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
    // Delete old events. We do not call clear() here since we want to
    // preserve any existing response function.
    if (m_events != NULL) delete m_events;

    // Allocate events
    GCTAEventCube* events = new GCTAEventCube;
    m_events = events;

    // Load events into cube
    events->load(filename);

    // Copy energy boundaries and GTIs from event cube
    m_ebounds = ((GCTAEventCube*)m_events)->m_ebds;
    m_gti     = ((GCTAEventCube*)m_events)->m_gti;

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
 *
 * @todo We allocate void response and pointing instances so make sure they
 * exist (analysis methods depend on the existence of these members).
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
 * @brief Temporally integrate spatially & spectrally integrated Npred
 *        gradient kernel
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
