/***************************************************************************
 *                          GTime.cpp - Time class                         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
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
#include "GTools.hpp"

/* __ Constants __________________________________________________________ */
const double sec_in_day = 86400.0;                     // Seconds in TT day
//const double mjd_ref    = 51910.0;                     // MJD of MET=0
//const double jd_ref     = 2451910.5;                   // JD seconds of MET=0
const double mjd_ref    = 51910.0007428703703703703;   // MJD of Fermi MET=0
const double jd_ref     = 2451910.5007428703703703703; // JD seconds of Fermi MET=0

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
 * @brief Return time in JD (TT) (unit: days)
 ***************************************************************************/
double GTime::jd(void) const
{
    // Convert time from MET to JD
    double jd = m_time / sec_in_day + jd_ref;
    
    // Return JD
    return jd;
}


/***********************************************************************//**
 * @brief Return time in MJD (unit: days)
 ***************************************************************************/
double GTime::mjd(void) const
{
    // Convert time from MET to MJD
    double mjd = m_time / sec_in_day + mjd_ref;
    
    // Return MJD
    return mjd;
}


/***********************************************************************//**
 * @brief Return time in MET (unit: seconds)
 ***************************************************************************/
double GTime::met(void) const
{
    // Return time
    return m_time;
}


/***********************************************************************//**
 * @brief Set time in JD (TT) (unit: days)
 *
 * @param[in] time Time in JD (TT).
 ***************************************************************************/
void GTime::jd(const double& time)
{
    // Convert time from JD to MET
    m_time = (time - jd_ref) * sec_in_day;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set time in MJD (TT) (unit: days)
 *
 * @param[in] time Time in MJD (TT).
 ***************************************************************************/
void GTime::mjd(const double& time)
{
    // Convert time from MJD to MET
    m_time = (time - mjd_ref) * sec_in_day;
    
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
    m_time = time;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Print time
 ***************************************************************************/
std::string GTime::print(void) const
{
    // Initialise result string
    std::string result;

    // Append time
    result.append(str(met())+" s");

    // Return
    return result;
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
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] time Time.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GTime& time)
{
     // Write time in output stream
    os << time.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] time Time.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GTime& time)
{
    // Write time into logger
    log << time.print();

    // Return logger
    return log;
}
