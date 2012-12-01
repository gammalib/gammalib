/***************************************************************************
 *                          GTime.cpp - Time class                         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
 * @brief Time class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTime.hpp"
#include "GTools.hpp"
#include "GException.hpp"

/* __ Constants __________________________________________________________ */
const double sec_in_day = 86400.0;                     // Seconds in TT day
//const double mjd_ref    = 51910.0;                     // MJD of MET=0
//const double jd_ref     = 2451910.5;                   // JD seconds of MET=0
const double mjd_ref    = 51910.0007428703703703703;   // MJD of Fermi MET=0
const double jd_ref     = 2451910.5007428703703703703; // JD seconds of Fermi MET=0

/* __ Method name definitions ____________________________________________ */
#define G_TIME                 "GTime::time(double&, double&, std::string&, " \
                                                 "std::string&, std::string&"

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
 * @brief Time constructor
 *
 * @param[in] time Time value.
 * @param[in] mjdref Reference MJD (days).
 * @param[in] timeunit Time unit ("sec(s)", "day(s)").
 * @param[in] timesys Time system (ignored so far).
 * @param[in] timeref Time reference (ignored so far).
 *
 * Sets the time for arbitrary units and MJD reference date. The MJD
 * reference day is specified as floating point value.
 ***************************************************************************/
GTime::GTime(const double&      time,
             const double&      mjdref,
             const std::string& timeunit,
             const std::string& timesys,
             const std::string& timeref)
{
    // Initialise private members for clean destruction
    init_members();

    // Set time
    this->time(time, mjdref, timeunit, timesys, timeref);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Time constructor
 *
 * @param[in] time Time value.
 * @param[in] mjdrefi Integer part of reference MJD (days).
 * @param[in] mjdreff Fractional part of reference MJD (days).
 * @param[in] timeunit Time unit (sec, days).
 * @param[in] timesys Time system (TT).
 * @param[in] timeref Local time reference.
 *
 * Sets the time for arbitrary units and MJD reference date. The MJD
 * reference day is specified as integer part and floating point fraction.
 ***************************************************************************/
GTime::GTime(const double&      time,
             const int&         mjdrefi,
             const double&      mjdreff,
             const std::string& timeunit,
             const std::string& timesys,
             const std::string& timeref)
{
    // Initialise private members for clean destruction
    init_members();

    // Set time
    this->time(time, mjdrefi, mjdreff, timeunit, timesys, timeref);

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
 * @brief Clear time
 ***************************************************************************/
void GTime::clear(void)
{
    // Free members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone object
 ***************************************************************************/
GTime* GTime::clone(void) const
{
    // Clone this image
    return new GTime(*this);
}


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
    // Convert time to MJD
    double mjd = m_time / sec_in_day + m_mjdref;
    
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
 * @brief Set time
 *
 * @param[in] time Time value.
 * @param[in] mjdref Reference MJD (days).
 * @param[in] timeunit Time unit ("s", "d", "sec(s)", "day(s)").
 * @param[in] timesys Time system (ignored so far).
 * @param[in] timeref Time reference (ignored so far).
 *
 * @exception GException::time_invalid_unit
 *            Invalid valid for "timeunit" encountered.
 *
 * Sets the time for arbitrary units and MJD reference date. The MJD
 * reference day is specified as floating point value.
 *
 * @todo Implement interpretation of "timesys" and "timeref" parameters.
 ***************************************************************************/
void GTime::time(const double&      time,
                 const double&      mjdref,
                 const std::string& timeunit,
                 const std::string& timesys,
                 const std::string& timeref)
{
    // Convert time to seconds
    double time_in_seconds;
    std::string u_timeunit = toupper(timeunit);
    if (u_timeunit == "D" || u_timeunit == "DAY" || u_timeunit == "DAYS") {
        time_in_seconds = time * sec_in_day;
    }
    else if (u_timeunit == "S" || u_timeunit == "SEC" || u_timeunit == "SECS") {
        time_in_seconds = time;
    }
    else {
        throw GException::time_invalid_unit(G_TIME, timeunit,
              "Valid units are \"d(ay(s))\" and \"s(ec(s))\"");
    }

    // Store time and reference
    m_time   = time_in_seconds;
    m_mjdref = mjdref;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set time
 *
 * @param[in] time Time value.
 * @param[in] mjdrefi Integer part of reference MJD (days).
 * @param[in] mjdreff Fractional part of reference MJD (days).
 * @param[in] timeunit Time unit (sec, days).
 * @param[in] timesys Time system (TT).
 * @param[in] timeref Local time reference.
 *
 * Sets the time for arbitrary units and MJD reference date. The MJD
 * reference day is specified as integer part and floating point fraction.
 ***************************************************************************/
void GTime::time(const double&      time,
                 const int&         mjdrefi,
                 const double&      mjdreff,
                 const std::string& timeunit,
                 const std::string& timesys,
                 const std::string& timeref)
{
    // Compute reference MJD
    double mjdref = double(mjdrefi) + mjdreff;
    
    // Set time
    this->time(time, mjdref, timeunit, timesys, timeref);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get time in seconds
 *
 * Returns time in seconds in the native MJD reference.
 ***************************************************************************/
double GTime::time(void) const
{
    // Return time
    return m_time;
}


/***********************************************************************//**
 * @brief Get time in seconds for a specific MJD
 *
 * @param[in] mjdref Reference MJD (days).
 *
 * Returns time in seconds in the specified MJD reference.
 ***************************************************************************/
double GTime::time(const double& mjdref) const
{
    // Retrieve time in seconds
    double time = m_time;
    
    // Add offset due to MJD differences
    if (m_mjdref != mjdref) {

        // Compute time offset
        double offset = (m_mjdref - mjdref) * sec_in_day;
        
        // Add offset
        time += offset;

    }

    // Return time
    return time;
}


/***********************************************************************//**
 * @brief Get time in seconds for a specific MJD
 *
 * @param[in] mjdrefi Integer part of reference MJD (days).
 * @param[in] mjdreff Fractional part of reference MJD (days).
 *
 * Returns time in seconds in the specified MJD reference.
 ***************************************************************************/
double GTime::time(const int& mjdrefi, const double& mjdreff) const
{
    // Compute reference MJD
    double mjdref = double(mjdrefi) + mjdreff;
    
    // Return time
    return (time(mjdref));
}


/***********************************************************************//**
 * @brief Get MJD reference in days
 *
 * Returns MJD reference in days.
 ***************************************************************************/
double GTime::mjdref(void) const
{
    // Return MJD reference
    return m_mjdref;
}


/***********************************************************************//**
 * @brief Print time
 *
 * Prints time in seconds in the native MJD reference.
 ***************************************************************************/
std::string GTime::print(void) const
{
    // Initialise result string
    std::string result;

    // Append time
    result.append(str(time())+" s");

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
 *
 * Initialises time to 0 and MJD reference to 2001-01-01 00:00:00.000 UTC.
 ***************************************************************************/
void GTime::init_members(void)
{
    // Initialise members
    m_time   = 0.0;
    m_mjdref = mjd_ref;
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] time Time.
 *
 * Copies the time and the MJD reference.
 ***************************************************************************/
void GTime::copy_members(const GTime& time)
{
    // Copy time
    m_time   = time.m_time;
    m_mjdref = time.m_mjdref;
    
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
