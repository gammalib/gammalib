/***************************************************************************
 *        GMWLObservation.cpp  -  Multi-wavelength observation class       *
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
 * @file GMWLObservation.cpp
 * @brief Multi-wavelength observation class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
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
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename File name.
 *
 * Creates instance from file.
 ***************************************************************************/
GMWLObservation::GMWLObservation(const std::string& filename) : GObservation()
{
    // Initialise members
    init_members();

    // Load observation
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename File name.
 * @param[in] extname FITS file extension name.
 *
 * Creates instance from file.
 ***************************************************************************/
GMWLObservation::GMWLObservation(const std::string& filename,
                                 const std::string& extname) : GObservation()
{
    // Initialise members
    init_members();

    // Load observation
    load(filename, extname);

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
    // Initialise members
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
 * @brief Returns pointer to response function (dummy)
 ***************************************************************************/
GMWLResponse* GMWLObservation::response(void) const
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
    name("Multi-wavelength observation");
    instrument(spec->instrument());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load observation
 *
 * @param[in] filename File name.
 * @param[in] extname FITS extension name.
 ***************************************************************************/
void GMWLObservation::load(const std::string& filename,
                           const std::string& extname)
{
    // Clear observation
    clear();

    // Allocate spectrum
    GMWLSpectrum* spec = new GMWLSpectrum;
    m_events = spec;

    // Load spectrum
    spec->load(filename, extname);

    // Set attributes
    name("Multi-wavelength observation");
    instrument(spec->instrument());

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
    result.append("=== GMWLObservation ===");
    result.append("\n"+parformat("Name")+name());

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
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
