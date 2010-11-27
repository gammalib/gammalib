/***************************************************************************
 *          GMWLDatum.cpp  -  Multi-wavelength spectral point class        *
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
 * @file GMWLDatum.cpp
 * @brief GMWLDatum class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include "GLog.hpp"
#include "GTools.hpp"
#include "GMWLDatum.hpp"

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
 * Creates instance of an undefined spectral point.
 ***************************************************************************/
GMWLDatum::GMWLDatum(void) : GEventBin()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] datum Spectral point.
 *
 * Creates instance by copying an existing spectral point.
 ***************************************************************************/
GMWLDatum::GMWLDatum(const GMWLDatum& datum) : GEventBin(datum)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(datum);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 *
 * Destroy instance.
 ***************************************************************************/
GMWLDatum::~GMWLDatum(void)
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
 * @param[in] datum Spectral point.
 *
 * Copies spectral point into the instance.
 ***************************************************************************/
GMWLDatum& GMWLDatum::operator= (const GMWLDatum& datum)
{
    // Execute only if object is not identical
    if (this != &datum) {

        // Copy base class members
        this->GEventBin::operator=(datum);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(datum);

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
 *
 * This method properly resets the instance to an initial state.
 ***************************************************************************/
void GMWLDatum::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GEventBin::free_members();
    this->GEvent::free_members();

    // Initialise members
    this->GEvent::init_members();
    this->GEventBin::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GMWLDatum* GMWLDatum::clone(void) const
{
    return new GMWLDatum(*this);
}


/***********************************************************************//**
 * @brief Print spectral point information
 ***************************************************************************/
std::string GMWLDatum::print(void) const
{
    // Initialise result string
    std::string result;

    // Append spectral point
    result.append(m_eng.print());
    if (m_eng_err.MeV() > 0.0)
        result.append(" +/- "+m_eng_err.print());
    result.append(": "+str(m_flux));
    if (m_flux_err > 0.0)
        result.append(" +/- "+str(m_flux_err));

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
void GMWLDatum::init_members(void)
{
    // Initialise members
    m_dir.clear();
    m_time.clear();
    m_eng.clear();
    m_eng_err.clear();
    m_flux     = 0.0;
    m_flux_err = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] datum Instance to be copied.
 ***************************************************************************/
void GMWLDatum::copy_members(const GMWLDatum& datum)
{
    // Copy members
    m_dir      = datum.m_dir;
    m_time     = datum.m_time;
    m_eng      = datum.m_eng;
    m_eng_err  = datum.m_eng_err;
    m_flux     = datum.m_flux;
    m_flux_err = datum.m_flux_err;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 *
 * This class does not allocate any memory but simply holds pointers. Hence
 * nothing has to be deallocated.
 ***************************************************************************/
void GMWLDatum::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] datum Spectral point.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GMWLDatum& datum)
{
     // Write spectral point in output stream
    os << datum.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] datum Spectral point.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GMWLDatum& datum)
{
    // Write spectral point into logger
    log << datum.print();

    // Return logger
    return log;
}
