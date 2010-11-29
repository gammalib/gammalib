/***************************************************************************
 *        GMWLObservation.cpp  -  Multi-wavelength observation class       *
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
 * @file GMWLObservation.cpp
 * @brief GMWLObservation class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include "GLog.hpp"
#include "GTools.hpp"
#include "GMWLObservation.hpp"
#include "GMWLSpectrum.hpp"

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
 * @brief Void constructor
 *
 * Creates instance of an undefined observation.
 ***************************************************************************/
GMWLObservation::GMWLObservation(void) : GObservation()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename Filename.
 *
 * Creates instance from file.
 ***************************************************************************/
GMWLObservation::GMWLObservation(const std::string& filename) : GObservation()
{
    // Initialise class members for clean destruction
    init_members();

    // Load observation
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] obs Observation.
 *
 * Creates instance by copying an existing observation.
 ***************************************************************************/
GMWLObservation::GMWLObservation(const GMWLObservation& obs) : GObservation(obs)
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
 *
 * Destroy instance.
 ***************************************************************************/
GMWLObservation::~GMWLObservation(void)
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
 * Copies observation into the instance.
 ***************************************************************************/
GMWLObservation& GMWLObservation::operator= (const GMWLObservation& obs)
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
 *
 * This method properly resets the instance to an initial state.
 ***************************************************************************/
void GMWLObservation::clear(void)
{
    // Free class members (base and derived classes, derived class first)
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
GMWLObservation* GMWLObservation::clone(void) const
{
    return new GMWLObservation(*this);
}


/***********************************************************************//**
 * @brief Set response function (dummy)
 *
 * @param[in] rspname Name of response function.
 * @param[in] caldb Optional name of calibration database.
 *
 * This method does nothing.
 ***************************************************************************/
void GMWLObservation::response(const std::string& rspname, std::string caldb)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns pointer to response function (dummy)
 *
 * @param[in] time Time.
 ***************************************************************************/
GMWLResponse* GMWLObservation::response(const GTime& time) const
{
    // Return response pointer
    return m_response;
}


/***********************************************************************//**
 * @brief Returns pointer to pointing (dummy)
 *
 * @param[in] time Time.
 ***************************************************************************/
GMWLPointing* GMWLObservation::pointing(const GTime& time) const
{
    // Return pointing pointer
    return m_pointing;
}


/***********************************************************************//**
 * @brief Load observation
 *
 * @param[in] filename File name.
 * @param[in] extname FITS extension name.
 *
 * @todo Extract energy boundaries from spectrum.
 *
 * @todo Set good time intervals covering observation interval.
 ***************************************************************************/
void GMWLObservation::load(const std::string& filename)
{
    // Clear observation
    clear();

    // Allocate spectrum
    GMWLSpectrum* spec = new GMWLSpectrum;
    m_events = spec;

    // Load spectrum
    spec->load(filename);

    // Set attributes
    obsname("Multi-wavelength observation");
    instrument(spec->instrument());

    // Extract energy boundaries from spectrum

    // Set good time intervals

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print multi-wavelength information
 ***************************************************************************/
std::string GMWLObservation::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GMWLObservation ===\n");
    result.append(parformat("Name")+obsname()+"\n");
    result.append(parformat("Instrument")+instrument()+"\n");
    result.append(parformat("Time range")+"\n");
    result.append(parformat("Energy range"));

    // Append events
    if (m_events != NULL)
        result.append("\n"+((GMWLSpectrum*)m_events)->print());

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GMWLObservation::init_members(void)
{
    // Initialise members
    m_instrument.clear();
    m_response = new GMWLResponse;
    m_pointing = new GMWLPointing;

    // Overwrite base class statistics
    m_statistics = "Gaussian";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] obs Observation to be copied
 ***************************************************************************/
void GMWLObservation::copy_members(const GMWLObservation& obs)
{
    // Copy members
    m_instrument = obs.m_instrument;
    if (obs.m_response != NULL) m_response = obs.m_response->clone();
    if (obs.m_pointing != NULL) m_pointing = obs.m_pointing->clone();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GMWLObservation::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] obs Observation.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GMWLObservation& obs)
{
     // Write observation in output stream
    os << obs.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] obs Observation.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GMWLObservation& obs)
{
    // Write observation into logger
    log << obs.print();

    // Return logger
    return log;
}
