/***************************************************************************
 *               GLATObservation.cpp  -  LAT Observation class             *
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
 * @file GLATObservation.cpp
 * @brief GLATObservationclass implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include "GException.hpp"
#include "GLATObservation.hpp"
#include "GLATEventList.hpp"
#include "GLATEventCube.hpp"
#include "GLATRoi.hpp"
#include "GFits.hpp"

/* __ Method name definitions ____________________________________________ */

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
GLATObservation::GLATObservation(void) : GObservation()
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
GLATObservation::GLATObservation(const GLATObservation& obs) : GObservation(obs)
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
GLATObservation::~GLATObservation(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Operators                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] obs Observation which should be assigned.
 ***************************************************************************/
GLATObservation& GLATObservation::operator= (const GLATObservation& obs)
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
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clone object
***************************************************************************/
GLATObservation* GLATObservation::clone(void) const
{
    return new GLATObservation(*this);
}


/***********************************************************************//**
 * @brief Load LAT response
 *
 * @param[in] irfname Name of instrument response function.
 * @param[in] caldb Optional path to calibration database.
 *
 * @todo Response not yet loaded.
 ***************************************************************************/
void GLATObservation::response(const std::string& irfname, std::string caldb)
{
    // Delete old response function
    if (m_response != NULL) delete m_response;

    // Allocate new LAT response function
    m_response = new GLATResponse;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns pointer to LAT response function
 *
 * @param[in] time Time.
 *
 * Returns pointer to response function for a given time. As the response is
 * supposed not to vary during an observation, the time argument needs not to
 * be considered.
 ***************************************************************************/
GResponse* GLATObservation::response(const GTime& time) const
{
    // Return response pointer
    return m_response;
}


/***********************************************************************//**
 * @brief Returns pointer to LAT pointing direction
 *
 * @param[in] time Time.
 *
 * Returns pointer to pointing direction for a given time. As the pointing
 * direction is supposed not to vary during an observation, the time argument
 * needs not to be considered.
 ***************************************************************************/
GPointing* GLATObservation::pointing(const GTime& time) const
{
    // Return response pointer
    return m_pointing;
}


/***********************************************************************//**
 * @brief Returns instrument name
 ***************************************************************************/
std::string GLATObservation::instrument(void) const
{
    // Return instument name
    return ("LAT");
}


/***********************************************************************//**
 * @brief Load data for unbinned analysis
 *
 * @param[in] ft1name FT1 FITS filename.
 * @param[in] ft2name FT2 FITS filename.
 * @param[in] ltcube_name Lifetime cube FITS filename
 *
 * @todo So far nothing is done with the ft2 file and the ltcube file.
 *       Loading of the relevant information needs to be implemented.
 * @todo Implement proper GTI loading method that provides correct time
 *       conversion
 ***************************************************************************/
void GLATObservation::load_unbinned(const std::string& ft1name,
                                    const std::string& ft2name,
                                    const std::string& ltcube_name)
{
    // Free and initialise base class members
    this->GObservation::free_members();
    this->GObservation::init_members();

    // Free and initialise class members
    free_members();
    init_members();

    // Allocate LAT events
    GLATEventList* events = new GLATEventList;
    m_events = events;

    // Load LAT events from FT1 file
    events->load(ft1name);

    // Load GTIs from FT1 file
    m_gti.load(ft1name);

    // Link observations to events. This has to be done after loading since
    // loading initialises the GLATEventList object, hence resets the pointer
    // to the observation.
    //events->obs(this);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load data for binned analysis
 *
 * @param[in] cntmap_name Counts map or Source map FITS filename
 * @param[in] expmap_name Binned explosure map FITS filename
 * @param[in] ltcube_name Lifetime cube FITS filename
 *
 * @todo So far nothing is done with the expmap and the ltcube files.
 *       Approriate loading needs to be implemented.
 * @todo Implement proper GTI loading method that provides correct time
 *       conversion
 ***************************************************************************/
void GLATObservation::load_binned(const std::string& cntmap_name,
                                  const std::string& expmap_name,
                                  const std::string& ltcube_name)
{
    // Free and initialise base class members
    this->GObservation::free_members();
    this->GObservation::init_members();

    // Free and initialise class members
    free_members();
    init_members();

    // Allocate LAT events
    GLATEventCube* events = new GLATEventCube;
    m_events = events;

    // Load LAT events from counts map file
    events->load(cntmap_name);

    // Copy over energy boundaries from events cube
    m_ebounds = events->m_ebds;

    // Load GTIs from counts map file
    m_gti.load(cntmap_name);

    // Set mean time
    events->m_time = 0.5 * (m_gti.tstart() + m_gti.tstop());

    // Link observations to events. This has to be done after loading since
    // loading initialises the GLATEventList object, hence resets the pointer
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
void GLATObservation::init_members(void)
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
void GLATObservation::copy_members(const GLATObservation& obs)
{
    // Copy members
    if (obs.m_response != NULL) m_response = obs.m_response->clone();
    if (obs.m_pointing != NULL) m_pointing = obs.m_pointing->clone();

    // Update the back pointer to link observation the actual observation
    // to the events. This has to be done here since the events that were
    // copied do not yet know to which observation they belong.
    /*
    if (m_events->islist())
        ((GLATEventList*)m_events)->obs(this);
    else
        ((GLATEventCube*)m_events)->obs(this);
    */

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATObservation::free_members(void)
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
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put LAT observation in output stream
 *
 * @param[in] os Output stream into which the data will be dumped
 * @param[in] obs Observation to be dumped
 *
 * @todo Implement GLATResponse operator<<
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GLATObservation& obs)
{
    // Put observation in stream
    os.precision(3);
    os << "=== GLATObservation ===" << std::endl;
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
        os << " Region of interest ........: " << *((GLATRoi*)obs.m_roi) << std::endl;

    // Add GTIs to stream
    os << obs.m_gti << std::endl;

    // Add response to stream if it exists
    //if (obs.m_response != NULL)
    //    os << *(obs.m_response) << std::endl;

    // Add events to stream
    if (obs.m_events != NULL) {
        if (obs.m_events->islist())
            os << *((GLATEventList*)obs.m_events) << std::endl;
        else
            os << *((GLATEventCube*)obs.m_events) << std::endl;
    }

    // Return output stream
    return os;
}
