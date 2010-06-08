/***************************************************************************
 *                          GTime.cpp - Time class                         *
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
 * @file GTime.cpp
 * @brief Time value class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTime.hpp"

/* __ Constants __________________________________________________________ */
const double sec_in_day = 86400.0;  // Seconds in day
const double met_ref    = 51910.0;  // MJD of MET=0 

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
GTime::GTime(void)
{
    // Initialise private members for clean destruction
    init_members();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] time Object from which the instance should be built.
 ***************************************************************************/
GTime::GTime(const GTime& time)
{ 
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(time);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GTime::~GTime(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] time Object which should be assigned.
 ***************************************************************************/
GTime& GTime::operator= (const GTime& time)
{ 
    // Execute only if object is not identical
    if (this != &time) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(time);

    } // endif: object was not identical
  
    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return time in MJD (unit: days)
 ***************************************************************************/
double GTime::mjd(void) const
{
    // Return time
    return m_time;
}


/***********************************************************************//**
 * @brief Return time in MET (unit: seconds)
 ***************************************************************************/
double GTime::met(void) const
{
    // Convert time from MJD to MET
    double met = (m_time - met_ref) * sec_in_day;
    
    // Return MET
    return met;
}


/***********************************************************************//**
 * @brief Set time in MJD (unit: days)
 *
 * @param[in] time Time in MJD.
 ***************************************************************************/
void GTime::mjd(const double& time)
{
    // Set time
    m_time = time;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set time in MET (unit: seconds)
 *
 * @param[in] time Time in MET.
 ***************************************************************************/
void GTime::met(const double& time)
{
    // Set time
    m_time = time / sec_in_day + met_ref;
    
    // Return
    return;
}
 

/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GTime::init_members(void)
{
    // Initialise members
    m_time = 0.0;
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] time Object from which members which should be copied.
 ***************************************************************************/
void GTime::copy_members(const GTime& time)
{
    // Copy time
    m_time = time.m_time;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GTime::free_members(void)
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
 * @brief Put object in output stream
 *
 * @param[in] os Output stream into which the model will be dumped
 * @param[in] time Object to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GTime& time)
{
    // Put object in stream
    os << time.mjd();

    // Return output stream
    return os;
}
