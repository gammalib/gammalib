/***************************************************************************
 *           GObservation.cpp  -  Observation abstract base class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2009 by Jurgen Knodlseder                           *
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
    // Free members
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
    m_tstart   = 0.0;
    m_tstop    = 0.0;
    m_emin     = 0.0;
    m_emax     = 0.0;
    m_events   = NULL;
    m_response = NULL;
    m_gti      = GGti();
    m_models   = GModels();

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
    m_models     = obs.m_models;

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
 =                           GObservation friends                          =
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
            obs.m_tstart << " - " << obs.m_tstop << std::endl;
    os << " Energy range ..............: " << std::fixed << 
            obs.m_emin << " - " << obs.m_emax << " MeV" << std::endl;

    // Add event list to stream
    if (obs.m_events != NULL)
        os << *(obs.m_events) << std::endl;

    // Add GTIs to stream
    os << obs.m_gti;
    
    // Add models to stream
    os << obs.m_models;
        
    // Return output stream
    return os;
}


/*==========================================================================
 =                                                                         =
 =                   Other functions used by GObservation                  =
 =                                                                         =
 ==========================================================================*/
