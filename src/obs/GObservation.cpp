/***************************************************************************
 *           GObservation.cpp  -  Observation abstract base class          *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
/**
 * @file GObservation.cpp
 * @brief GObservation abstract base class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GException.hpp"
#include "GObservation.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_RELEASE      "GObservation::release(void)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                    GObservation constructors/destructors                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GObservation::GObservation()
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
GObservation::~GObservation()
{
	// Do not destruct if object is linked in somewhere
	if (m_num_links > 0)
		m_free_pending = 1;
	
    // Free members if not linked
	else
		free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                          GObservation operators                         =
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
 =                       GObservation public methods                       =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Link observation
 ***************************************************************************/
void GObservation::link(void)
{
	// Increment number of links
	m_num_links++;
	
    // Return
    return;
}


/***********************************************************************//**
 * @brief Release observation
 ***************************************************************************/
void GObservation::release(void)
{
	// Decrement number of links
	if (m_num_links > 0)
		m_num_links--;
	else
		throw GException::invalid_release(G_RELEASE, m_num_links);
	
	// If there are no more links and destruction is pending then free
	// object now
	if (m_num_links == 0 && m_free_pending == 1)
		free_members();
	
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                       GObservation private methods                      =
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
    m_tstart       = 0.0;
    m_tstop        = 0.0;
    m_events       = NULL;
    m_gti          = NULL;
    m_response     = NULL;
	m_num_links    = 0;
	m_free_pending = 0;

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
    m_obsname      = obs.m_obsname;
    m_instrument   = obs.m_instrument;
    m_tstart       = obs.m_tstart;
    m_tstop        = obs.m_tstop;
	m_num_links    = 0;           // Do not copy over any link
	m_free_pending = 0;           // Do not copy over destruction scheduling

    // Clone members
	m_events   = obs.m_events->clone();
	m_gti      = obs.m_gti->clone();
	m_response = obs.m_response->clone();

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
	if (m_gti      != NULL) delete m_gti;
	if (m_response != NULL) delete m_response;

    // Signal free pointers
    m_events   = NULL;
    m_gti      = NULL;
    m_response = NULL;
	
	// Reset number of links and destruction scheduling
	m_num_links    = 0;
	m_free_pending = 0;

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                           GObservation friends                          =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                   Other functions used by GObservation                  =
 =                                                                         =
 ==========================================================================*/
