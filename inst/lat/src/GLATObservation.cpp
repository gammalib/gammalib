/***************************************************************************
 *               GLATObservation.cpp  -  LAT Observation class             *
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
 * @file GLATObservation.cpp
 * @brief LAT Observation class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GLATObservation.hpp"
#include "GLATEventList.hpp"
#include "GLATEventCube.hpp"
#include "GLATRoi.hpp"
#include "GFits.hpp"
#include "GTools.hpp"
#include "GEnergy.hpp"

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
 * @brief Void constructor
 ***************************************************************************/
GLATObservation::GLATObservation(void) : GObservation()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] obs LAT observation.
 ***************************************************************************/
GLATObservation::GLATObservation(const GLATObservation& obs) : GObservation(obs)
{
    // Initialise members
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
 * @param[in] obs LAT observation.
 ***************************************************************************/
GLATObservation& GLATObservation::operator= (const GLATObservation& obs)
{
    // Execute only if object is not identical
    if (this != &obs) {

        // Copy base class members
        this->GObservation::operator=(obs);

        // Free members
        free_members();

        // Initialise members
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
 * @brief Clear instance
 ***************************************************************************/
void GLATObservation::clear(void)
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

    // Set calibration database
    m_response->caldb(caldb);

    // Load instrument response function
    m_response->load(irfname);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns pointer to LAT response function
 ***************************************************************************/
GLATResponse* GLATObservation::response(void) const
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
GLATPointing* GLATObservation::pointing(const GTime& time) const
{
    // Return response pointer
    return m_pointing;
}


/***********************************************************************//**
 * @brief Returns pointer to LAT livetime cube
 ***************************************************************************/
GLATLtCube* GLATObservation::ltcube(void) const
{
    // Return livetime cube pointer
    return m_ltcube;
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
 * @brief Print LAT observation information
 ***************************************************************************/
std::string GLATObservation::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GLATObservation ===\n");
    result.append(parformat("Name")+name()+"\n");
    result.append(parformat("Instrument")+instrument()+"\n");
    result.append(parformat("Statistics")+statistics()+"\n");

    // Append response
    if (m_response != NULL)
        result.append("\n"+m_response->print());
    else
        result.append("\n"+parformat("LAT response")+"undefined");

    // Append livetime cube
    if (m_ltcube != NULL)
        result.append("\n"+m_ltcube->print());
    else
        result.append("\n"+parformat("LAT livetime cube")+"undefined");

    // Append events
    if (m_events != NULL)
        result.append("\n"+m_events->print());

    // Return result
    return result;
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
 ***************************************************************************/
void GLATObservation::load_unbinned(const std::string& ft1name,
                                    const std::string& ft2name,
                                    const std::string& ltcube_name)
{
    // Delete any existing event container (do not call clear() as we do not
    // want to delete the response function)
    if (m_events != NULL) delete m_events;
    m_events = NULL;

    // Allocate event list
    GLATEventList* events = new GLATEventList;

    // Assign event list as the observation's event container
    m_events = events;

    // Open FITS file
    GFits file(ft1name);

    // Read event list
    events->read(file);

    // Read observation attributes from EVENTS extension
    //GFitsHDU* hdu = file.hdu("EVENTS");
    //read_attributes(hdu);

    // Close FITS file
    file.close();

    // Optionally allocate and load livetime cube
    if (ltcube_name.length() > 0) {
        m_ltcube = new GLATLtCube;
        m_ltcube->load(ltcube_name);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load data for binned analysis
 *
 * @param[in] cntmap_name Counts map or Source map FITS filename
 * @param[in] expmap_name Binned explosure map FITS filename
 * @param[in] ltcube_name Livetime cube FITS filename
 *
 * @todo So far nothing is done with the expmap file.
 *       Approriate loading needs to be implemented.
 ***************************************************************************/
void GLATObservation::load_binned(const std::string& cntmap_name,
                                  const std::string& expmap_name,
                                  const std::string& ltcube_name)
{
    // Delete old events and livetime cube.  We do not call clear() here
    // since we want to preserve any existing response function.
    if (m_events != NULL) delete m_events;
    if (m_ltcube != NULL) delete m_ltcube;
    m_events = NULL;
    m_ltcube = NULL;

    // Allocate event cube
    GLATEventCube* events = new GLATEventCube;

    // Assign event cube as the observation's event container
    m_events = events;

    // Load event list
    events->load(cntmap_name);

    // Optionally allocate and load livetime cube
    if (ltcube_name.length() > 0) {
        m_ltcube = new GLATLtCube;
        m_ltcube->load(ltcube_name);
    }

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
    m_ltcube   = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] obs LAT observation.
 ***************************************************************************/
void GLATObservation::copy_members(const GLATObservation& obs)
{
    // Copy members
    if (obs.m_response != NULL) m_response = obs.m_response->clone();
    if (obs.m_pointing != NULL) m_pointing = obs.m_pointing->clone();
    if (obs.m_ltcube   != NULL) m_ltcube   = obs.m_ltcube->clone();

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
    if (m_ltcube   != NULL) delete m_ltcube;

    // Mark memory as free
    m_response = NULL;
    m_pointing = NULL;
    m_ltcube   = NULL;

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/
