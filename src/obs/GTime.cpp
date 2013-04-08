/***************************************************************************
 *                          GTime.cpp - Time class                         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
const double mjd_ref = 55197.000766018518519;             //!< MJD of time=0
const double jd_ref  = mjd_ref + 2400000.5;               //!< JD of time=0

/* __ Method name definitions ____________________________________________ */
#define G_CONSTRUCT                     "GTime::GTime(double&, std::string&)"
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
    // Initialise private members
    init_members();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] time Time.
 ***************************************************************************/
GTime::GTime(const GTime& time)
{ 
    // Initialise private members
    init_members();

    // Copy members
    copy_members(time);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Time constructor
 *
 * @param[in] time Time value (seconds or days).
 * @param[in] unit Time unit (default: "secs").
 *
 * @exception GException::time_invalid_unit
 *            Invalid time unit specified.
 *
 * Constructs a GTime object by setting the time in the native reference
 * in units of seconds (default) or days.
 ***************************************************************************/
GTime::GTime(const double& time, const std::string& unit)
{
    // Initialise private members
    init_members();

    // Set time according to timeunit string
    std::string timeunit = tolower(unit);
    if (timeunit == "d" || timeunit == "day" || timeunit == "days") {
        days(time);
    }
    else if (timeunit == "s" || timeunit == "sec" || timeunit == "secs") {
        secs(time);
    }
    else {
        throw GException::time_invalid_unit(G_CONSTRUCT, unit,
              "Valid timeunit values are: \"d\", \"day\", \"days\","
              " \"s\", \"sec\" or \"secs\"");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Time constructor
 *
 * @param[in] time Time in given reference system.
 * @param[in] ref Reference system.
 *
 * Constructs a GTime object by setting the time to a value given in a
 * specific reference system.
 ***************************************************************************/
GTime::GTime(const double& time, const GTimeReference& ref)
{
    // Initialise private members
    init_members();

    // Set time
    set(time, ref);

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
 * @param[in] time Time.
 * @return Time.
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
 * @brief Clone time
 *
 * @return Pointer to deep copy of time.
 ***************************************************************************/
GTime* GTime::clone(void) const
{
    // Clone this image
    return new GTime(*this);
}


/***********************************************************************//**
 * @brief Return time in Julian Days (TT) (unit: days)
 *
 * @return Time in Julian Days (TT) [days].
 *
 * Returns the time in Julian Days (JD) in the Terrestrial Time (TT) system.
 ***************************************************************************/
double GTime::jd(void) const
{
    // Convert time from MET to JD
    double jd = m_time / sec_in_day + jd_ref;
    
    // Return JD
    return jd;
}


/***********************************************************************//**
 * @brief Return time in Modified Julian Days (TT) (unit: days)
 *
 * @return Time in Modified Julian Days (TT) [days].
 *
 * Returns the time in Modified Julian Days (MJD) in the Terrestrial Time
 * (TT) system.
 ***************************************************************************/
double GTime::mjd(void) const
{
    // Convert time to MJD
    double mjd = m_time / sec_in_day + mjd_ref;
    
    // Return MJD
    return mjd;
}


/***********************************************************************//**
 * @brief Return time in native reference (TT) (unit: seconds)
 *
 * @return Time in native reference [seconds].
 ***************************************************************************/
double GTime::secs(void) const
{
    // Return time
    return m_time;
}


/***********************************************************************//**
 * @brief Return time in native reference (TT) (unit: days)
 *
 * @return Time in native reference [days].
 ***************************************************************************/
double GTime::days(void) const
{
    // Return time
    return m_time / sec_in_day;
}


/***********************************************************************//**
 * @brief Return time in specified reference
 *
 * @return Time in specified reference.
 *
 * Convert the time from the native reference system into the specified
 * reference system.
 *
 * @todo Implement TT-UTC conversion if required. This requires
 *       implementation of leap seconds.
 ***************************************************************************/
double GTime::convert(const GTimeReference& ref) const
{
    // Retrieve time in native reference (seconds)
    double time = m_time;
    
    // Compute time offset in seconds
    double offset = (mjd_ref - ref.mjdref()) * sec_in_day;
        
    // Add time offset in seconds
    time += offset;

    // Convert to specified time unit
    double to_unit = ref.unitseconds();
    if (to_unit != 1.0) {
        time /= to_unit;
    }

    // Return time
    return time;
}


/***********************************************************************//**
 * @brief Set time in Julian Days (TT) (unit: days)
 *
 * @param[in] time Time in Julian Days (TT) [days].
 ***************************************************************************/
void GTime::jd(const double& time)
{
    // Convert time from JD to native (seconds)
    m_time = (time - jd_ref) * sec_in_day;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set time in Modified Julian Days (TT) (unit: days)
 *
 * @param[in] time Time in Modified Julian Days (TT) [days].
 ***************************************************************************/
void GTime::mjd(const double& time)
{
    // Convert time from MJD to native (seconds)
    m_time = (time - mjd_ref) * sec_in_day;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set time in native reference in seconds (TT)
 *
 * @param[in] seconds Time (TT) [seconds].
 ***************************************************************************/
void GTime::secs(const double& seconds)
{
    // Set time
    m_time = seconds;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set time in native reference in days (TT)
 *
 * @param[in] days Time (TT) [days].
 ***************************************************************************/
void GTime::days(const double& days)
{
    // Set time
    m_time = days * sec_in_day;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set time given in specified reference
 *
 * @param[in] time Time in given reference system.
 * @param[in] ref Reference system.
 *
 * Set the time to a value given in a specific reference system.
 *
 * @todo Implement TT-UTC conversion if required. This requires
 *       implementation of leap seconds.
 ***************************************************************************/
void GTime::set(const double& time, const GTimeReference& ref)
{
    // Convert time to specified time unit
    m_time = time * ref.unitseconds();
    
    // Compute time offset in seconds
    double offset = (mjd_ref - ref.mjdref()) * sec_in_day;
        
    // Subtract time offset in seconds
    m_time -= offset;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns native time reference
 *
 * @return Native time reference.
 *
 * Returns the native GammaLib time reference. The GammaLib native time
 * reference (i.e. time=0) is defined as January 1, 2010, 00:00:00 (TT).
 * The time system is Terrestrial Time (TT). Time is stored in seconds.
 ***************************************************************************/
GTimeReference GTime::reference(void) const
{
    // Allocate native time reference
    GTimeReference reference(mjd_ref, "s", "TT", "LOCAL");

    // Return reference
    return reference;
}


/***********************************************************************//**
 * @brief Print time
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing time in seconds in native reference.
 *
 * Prints time in seconds in the native reference.
 ***************************************************************************/
std::string GTime::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append time
        result.append(str(m_time)+" s (TT)");

    } // endif: chatter was not silent

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
 * @param[in] time Time.
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
