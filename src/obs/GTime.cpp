/***************************************************************************
 *                          GTime.cpp - Time class                         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2017 by Juergen Knoedlseder                         *
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
#include <ctime>
#include <cstring>      // std::memcpy
#include <cstdio>
#include "GTools.hpp"
#include "GMath.hpp"
#include "GException.hpp"
#include "GTime.hpp"
#include "GTimeReference.hpp"

/* __ Constants __________________________________________________________ */
const double mjd_ref = 55197.000766018518519;             //!< MJD of time=0
const double jd_ref  = mjd_ref + 2400000.5;               //!< JD of time=0

/* __ Method name definitions ____________________________________________ */
#define G_CONSTRUCT                     "GTime::GTime(double&, std::string&)"
#define G_SECS_GET                                "GTime::secs(std::string&)"
#define G_SECS_SET                       "GTime::secs(double&, std::string&)"
#define G_UTC                                      "GTime::utc(std::string&)"

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
 * @param[in] time Time value in TT (seconds or days).
 * @param[in] unit Time unit string.
 *
 * @exception GException::time_invalid_unit
 *            Invalid time unit specified.
 *
 * Constructs a GTime object by setting the time in the native reference
 * in the TT time system in units of seconds (default) or days.
 ***************************************************************************/
GTime::GTime(const double& time, const std::string& unit)
{
    // Initialise private members
    init_members();

    // Set time according to timeunit string
    std::string timeunit = gammalib::tolower(unit);
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
 * @brief Time constructor
 *
 * @param[in] time Time string in given reference system.
 * @param[in] ref Reference system.
 *
 * Constructs a GTime object by setting the time to a string value. See the
 * set(std::string&, GTimeReference&) method for valid time strings.
 ***************************************************************************/
GTime::GTime(const std::string& time, const GTimeReference& ref)
{
    // Initialise private members
    init_members();

    // Set time
    set(time, ref);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Time constructor
 *
 * @param[in] time Time string.
 *
 * Constructs a GTime object by setting the time to a string value. See the
 * set(const std::string& time) method for valid time strings.
 ***************************************************************************/
GTime::GTime(const std::string& time)
{
    // Initialise private members
    init_members();

    // Set time
    set(time);

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
GTime& GTime::operator=(const GTime& time)
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
 * @brief Return time in Julian Days (TT)
 *
 * @return Time in Julian Days (TT) (days).
 *
 * Returns the time in Julian Days (JD) in the Terrestrial Time (TT) system.
 ***************************************************************************/
double GTime::jd(void) const
{
    // Convert time from MET to JD
    double jd = m_time * gammalib::sec2day + jd_ref;
    
    // Return Julian Days
    return jd;
}


/***********************************************************************//**
 * @brief Return time in Julian Days for various time system
 *
 * @param[in] timesys Time system (one of "TT", "TAI", "UTC")
 * @return Time in Julian Days (days).
 *
 * Returns the time in Julian Days (JD) for the specified time system.
 ***************************************************************************/
double GTime::jd(const std::string& timesys) const
{
    // Convert time to Julian Days
    double jd = secs(timesys) * gammalib::sec2day + jd_ref;
    
    // Return Julian Days
    return jd;
}


/***********************************************************************//**
 * @brief Return time in Modified Julian Days (TT)
 *
 * @return Time in Modified Julian Days (TT) (days).
 *
 * Returns the time in Modified Julian Days (MJD) in the Terrestrial Time
 * (TT) system.
 ***************************************************************************/
double GTime::mjd(void) const
{
    // Convert time to Modified Julian Days
    double mjd = m_time * gammalib::sec2day + mjd_ref;
    
    // Return Modified Julian Days
    return mjd;
}


/***********************************************************************//**
 * @brief Return time in Modified Julian Days for various time system
 *
 * @param[in] timesys Time system (one of "TT", "TAI", "UTC")
 * @return Time in Modified Julian Days (days).
 *
 * Returns the time in Modified Julian Days (JD) for the specified time
 * system.
 ***************************************************************************/
double GTime::mjd(const std::string& timesys) const
{
    // Convert time to Modified Julian Days
    double mjd = secs(timesys) * gammalib::sec2day + mjd_ref;
    
    // Return Modified Julian Days
    return mjd;
}


/***********************************************************************//**
 * @brief Return time in seconds in native reference for various time systems
 *
 * @param[in] timesys Time system (one of "TT", "TAI", "UTC")
 * @return Time (seconds).
 ***************************************************************************/
double GTime::secs(const std::string& timesys) const
{
    // Initialise time
    double time;

    // If system is TT then return the stored time ...
    if (timesys == "TT") {
        time = m_time;
    }
    
    // ... otherwise if the system if TAI then subtract an offset ...
    else if (timesys == "TAI") {
        time = m_time - gammalib::tai2tt;
    }

    // ... otherwise if the system is UTC then convert the time from TT to TAI
    // and subtract the leap seconds. The repeated calling of the leap second
    // method is a kluge to converge as the argument is given in the UTC system
    // but upon start we're in the TAI system.
    else if (timesys == "UTC") {
        time         = m_time - gammalib::tai2tt; // Time in TAI
        double mjd   = time * gammalib::sec2day + mjd_ref;
        double leaps = leap_seconds(mjd);
        leaps        = leap_seconds(mjd - leaps * gammalib::sec2day);
        leaps        = leap_seconds(mjd - leaps * gammalib::sec2day);
        time        -= leaps;
    }

    // ... otherwise throw an exception
    else {
        std::string msg = "Unknown time system \""+timesys+"\". Either specify "
                          "\"TT\", \"TAI\" or \"UTC\".";
        throw GException::invalid_argument(G_SECS_GET, msg);
    }

    // Return time
    return time;
}


/***********************************************************************//**
 * @brief Return time in days in native reference (TT)
 *
 * @return Time in native reference (days).
 ***************************************************************************/
double GTime::days(void) const
{
    // Return time
    return (m_time * gammalib::sec2day);
}


/***********************************************************************//**
 * @brief Return time in days in native reference for various time system
 *
 * @param[in] timesys Time system (one of "TT", "TAI", "UTC")
 * @return Time in native reference (days).
 ***************************************************************************/
double GTime::days(const std::string& timesys) const
{
    // Return time
    return (secs(timesys) * gammalib::sec2day);
}


/***********************************************************************//**
 * @brief Return time as string in UTC time system
 *
 * @return Time as string in UTC time system.
 *
 * Returns time in the format YYYY-MM-DDThh:mm:ss, where YYYY is a four-digit
 * year, MM a two-digit month, DD a two-digit day of month, hh two digits of
 * hour (0 through 23), mm two digits of minutes, and ss two digits of
 * second (ISO 8601 time standard).
 *
 * The method is only valid for dates from year 1972 on.
 ***************************************************************************/
std::string GTime::utc(void) const
{
    // Define number of days per month
    static int daymonth[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    // Get MJD in TAI time system
    double mjd = this->mjd() - gammalib::tai2tt * gammalib::sec2day;
    
    // Correct for leap seconds (repeat is a kluge to converge as the
    // argument is given in the UTC system but upon start we're in the TAI
    // system
    double leaps = leap_seconds(mjd) * gammalib::sec2day;
    leaps        = leap_seconds(mjd - leaps) * gammalib::sec2day;
    leaps        = leap_seconds(mjd - leaps) * gammalib::sec2day;
    mjd         -= leaps;

    // Split in day and fraction
    int    day      = (int)mjd;
    double fraction = mjd - (double)day;

    // Compute time in day. We add a margin of 0.5 to the seconds and
    // subtract it later to avoid rounding of 59 to 60.
    double second = fraction * gammalib::sec_in_day + 0.5;
    int    hour   = (int)second / 3600;
    second       -= hour * 3600.0;
    int minute    = (int)second / 60;
    second       -= minute * 60.0;
    if (hour > 23) {
        hour -= 24;
        day++;
    }
    second -= 0.5;
    if (second < 0.0) {
        second = 0.0;
    }

    // Compute year and day in the year
    int year = 1972;           // Set year to 1972
    day     -= 41317;          // Subtract MJD of 1-1-1972
    day++;                     // Day 0 is 1 January
    int days = days_in_year(year);
    while (day > days) {
        day -= days;
        year++;
        days = days_in_year(year);
    }

    // Adjust number of days per month for leap years
    if (is_leap_year(year)) {
        daymonth[1] = 29;
    }
    else {
        daymonth[1] = 28;
    }

    // Compute month and day in the month
    int month = 0;
    while (month < 12) {
        if (day <= daymonth[month]) {
            break;
        }
        day -= daymonth[month];
        month++;
    }
    month++;

    // Create string
    char utc[32];
    sprintf(utc, "%4.4d-%2.2d-%2.2dT%2.2d:%2.2d:%02.06f",
            year, month, day, hour, minute, second);

    // Return
    return (std::string(utc));
}


/***********************************************************************//**
 * @brief Return Greenwich mean sidereal time in hours in a day
 *
 * @return Greenwich mean sidereal time (hours).
 *
 * See http://aa.usno.navy.mil/faq/docs/GAST.php
 ***************************************************************************/
double GTime::gmst(void) const
{
    // Get days since 1 January 2000, 12h in Universal Time
    double d = days("UTC") + 3652.50076602;

    // Compute Greenwich mean sidereal time in hours
    double gmst = 18.697374558 + 24.06570982441908 * d;

    // Put into [0,24]
    gmst -= floor(gmst/24.0) * 24.0;

    // Return Greenwich mean sidereal time in hours
    return gmst;
}


/***********************************************************************//**
 * @brief Return Greenwich apparent sidereal time in hours in a day
 *
 * @return Greenwich apparent sidereal time (hours).
 *
 * See http://aa.usno.navy.mil/faq/docs/GAST.php
 ***************************************************************************/
double GTime::gast(void) const
{
    // Get days since 1 January 2000, 12h in Universal Time
    double d = days("UTC") + 3652.50076602;

    // Compute longitude of the ascending node of the Moon in degrees
    double Omega = 125.04 - 0.052954 * d;

    // Compute mean longitude of the Sun in degrees
    double L = 280.47 + 0.98565 * d;

    // Compute nutation in longitude in hours
    double DeltaPsi = -0.000319 * std::sin(Omega * gammalib::deg2rad) -
                       0.000024 * std::sin(2.0 * L * gammalib::deg2rad);

    // Compute the obliquity in degrees
    double epsilon = 23.4393 - 0.0000004 * d;

    // Compute equation of equinoxes
    double eqeq = DeltaPsi * std::cos(epsilon * gammalib::deg2rad);

    // Compute Greenwich apparent sidereal time in hours
    double gast = gmst() + eqeq;

    // Put into [0,24]
    gast -= floor(gast/24.0) * 24.0;

    // Return Greenwich apparent sidereal time in hours
    return gast;
}


/***********************************************************************//**
 * @brief Return local mean sidereal time in hours in a day
 *
 * @param[in] geolon Geographic longitude West of Greenwich (degrees).
 * @return Local mean sidereal time (hours).
 *
 * See http://aa.usno.navy.mil/faq/docs/GAST.php
 ***************************************************************************/
double GTime::lmst(const double& geolon) const
{
    // Compute local mean siderial time
    double lmst = gmst() - geolon/15.0;

    // Put into [0,24]
    lmst -= floor(lmst/24.0) * 24.0;

    // Return local mean sidereal time in hours
    return lmst;
}


/***********************************************************************//**
 * @brief Return local apparent sidereal time in hours in a day
 *
 * @param[in] geolon Geographic longitude West of Greenwich (degrees).
 * @return Local apparent sidereal time (hours).
 *
 * See http://aa.usno.navy.mil/faq/docs/GAST.php
 ***************************************************************************/
double GTime::last(const double& geolon) const
{
    // Compute local mean siderial time
    double last = gast() - geolon/15.0;

    // Put into [0,24]
    last -= floor(last/24.0) * 24.0;

    // Return local mean sidereal time in hours
    return last;
}


/***********************************************************************//**
 * @brief Return time in specified reference
 *
 * @return Time in specified reference.
 *
 * Convert the time from the native reference system into the specified
 * reference system.
 ***************************************************************************/
double GTime::convert(const GTimeReference& ref) const
{
    // Retrieve time in native reference (TT in seconds)
    double time = m_time;
    
    // Compute time offset in seconds
    double offset = (mjd_ref - ref.mjdref()) * gammalib::sec_in_day;
        
    // Add time offset in seconds
    time += offset;

    // Subtract leap seconds in case that time is requested in UTC time
    // system
    if (gammalib::toupper(ref.timesys()) == "UTC") {

        // Get MJD in TAI time system
        double mjd = this->mjd() - gammalib::tai2tt * gammalib::sec2day;
    
        // Get leap seconds (repeat is a kluge to converge as the
        // argument is given in the UTC system but upon start we're in the TAI
        // system
        double leaps = leap_seconds(mjd);
        leaps        = leap_seconds(mjd - leaps * gammalib::sec2day);
        leaps        = leap_seconds(mjd - leaps * gammalib::sec2day);

        // Subtract leap seconds and TAI offset
        time -= (leaps + gammalib::tai2tt);

    } // endif: UTC time system has been requested

    // ... otherwise, if time is requested in TAI system then subtract the
    // TAI offset
    else if (gammalib::toupper(ref.timesys()) == "TAI") {
        time -= gammalib::tai2tt;
    }

    // Convert to specified time unit
    double to_unit = ref.unitseconds();
    if (to_unit != 1.0) {
        time /= to_unit;
    }

    // Return time
    return time;
}


/***********************************************************************//**
 * @brief Set time in Julian Days in native reference (TT)
 *
 * @param[in] time Time in Julian Days (TT) (days).
 ***************************************************************************/
void GTime::jd(const double& time)
{
    // Convert time from Julian Days to seconds in native reference
    m_time = (time - jd_ref) * gammalib::sec_in_day;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set time in Julian Days in native reference for various time systems
 *
 * @param[in] time Time in Julian Days (days).
 * @param[in] timesys Time system (one of "TT", "TAI", "UTC")
 ***************************************************************************/
void GTime::jd(const double& time, const std::string& timesys)
{
    // Convert time from Julian Days to seconds in native reference
    double seconds = (time - jd_ref) * gammalib::sec_in_day;

    // Set time according to the specified time system
    secs(seconds, timesys);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set time in Modified Julian Days in native reference (TT)
 *
 * @param[in] time Time in Modified Julian Days (TT) (days).
 ***************************************************************************/
void GTime::mjd(const double& time)
{
    // Convert time from Modified Julian Days to native (seconds)
    m_time = (time - mjd_ref) * gammalib::sec_in_day;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set time in Modified Julian Days in native reference for various
 *        time systems
 *
 * @param[in] time Time in Modified Julian Days (days).
 * @param[in] timesys Time system (one of "TT", "TAI", "UTC")
 ***************************************************************************/
void GTime::mjd(const double& time, const std::string& timesys)
{
    // Convert time from Modified Julian Days to seconds in native reference
    double seconds = (time - mjd_ref) * gammalib::sec_in_day;

    // Set time according to the specified time system
    secs(seconds, timesys);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set time in seconds in native reference for various time systems
 *
 * @param[in] seconds Time in native reference (seconds).
 * @param[in] timesys Time system (one of "TT", "TAI", "UTC")
 ***************************************************************************/
void GTime::secs(const double& seconds, const std::string& timesys)
{
    // If system is TT then simply set time ...
    if (timesys == "TT") {
        m_time = seconds;
    }
    
    // ... otherwise if the system if TAI then add an offset ...
    else if (timesys == "TAI") {
        m_time = seconds + gammalib::tai2tt;
    }

    // ... otherwise if the system is UTC then convert the time from TT to TAI
    // and add the leap seconds ...
    else if (timesys == "UTC") {
        double mjd   = seconds * gammalib::sec2day + mjd_ref; // Time in UTC
        double leaps = leap_seconds(mjd);
        m_time       = seconds + gammalib::tai2tt + leaps;
    }

    // ... otherwise throw an exception
    else {
        std::string msg = "Unknown time system \""+timesys+"\". Either specify "
                          "\"TT\", \"TAI\" or \"UTC\".";
        throw GException::invalid_argument(G_SECS_SET, msg);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set time in days in native reference (TT)
 *
 * @param[in] days Time (TT) (days).
 ***************************************************************************/
void GTime::days(const double& days)
{
    // Set time
    m_time = days * gammalib::sec_in_day;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set time in days in native reference for various time systems
 *
 * @param[in] days Time (TT) (days).
 * @param[in] timesys Time system (one of "TT", "TAI", "UTC")
 ***************************************************************************/
void GTime::days(const double& days, const std::string& timesys)
{
    // Convert time from days to seconds in native reference
    double seconds = days * gammalib::sec_in_day;

    // Set time according to the specified time system
    secs(seconds, timesys);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set time as string in UTC time system
 *
 * @param[in] time Time string (UTC).
 *
 * @exception GException::invalid_argument
 *            Invalid time string specified.
 *
 * The time has to be given in the format YYYY-MM-DDThh:mm:ss.s, where
 * YYYY is a four-digit year, MM a two-digit month, DD a two-digit day of
 * month, hh two digits of hour (0 through 23), mm two digits of minutes,
 * ss two digits of second and s one or more digits representing a
 * decimal fraction of a second (ISO 8601 time standard).
 *
 * The method is only valid for dates from year 1972 on.
 ***************************************************************************/
void GTime::utc(const std::string& time)
{
    // Define number of days per month
    static int daymonth[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    // Analyse the UTC string
    long   year   = 0;
    int    month  = 0;
    long   day    = 0;
    long   hour   = 0;
    long   minute = 0;
    double second = 0.0;
    int    n      = sscanf(time.c_str(), "%ld-%d-%ldT%ld:%ld:%lg",
		                   &year, &month, &day, &hour, &minute, &second);

    // Adjust number of days per month for leap years
    if (is_leap_year(year)) {
        daymonth[1] = 29;
    }
    else {
        daymonth[1] = 28;
    }

    // Check time string
    if (n != 3 && n != 6) {
        std::string msg = "Invalid time string \""+time+"\" encountered. "
                          "Please specify the time in the format YYYY-MM-DD "
                          "or YYYY-MM-DDThh:mm:ss.s.";
        throw GException::invalid_argument(G_UTC, msg);
    }
    if (year < 1000 || year > 9999) {
        std::string msg = "Invalid year "+gammalib::str(year)+" specified. "
                          "The year needs to be a four-digit year.";
        throw GException::invalid_argument(G_UTC, msg);
    }
    if (month < 1 || month > 12) {
        std::string msg = "Invalid month "+gammalib::str(month)+" specified. "
                          "The month needs to be comprised between 01 and 12";
        throw GException::invalid_argument(G_UTC, msg);
    }
    if (day < 1 || day > daymonth[month-1]) {
        std::string msg = "Invalid day "+gammalib::str(day)+" specified. "
                          "The day for month "+gammalib::str(month)+" needs "
                          "to be comprised between 01 and "+
                          gammalib::str(daymonth[month-1])+".";
        throw GException::invalid_argument(G_UTC, msg);
    }
    if (hour < 0 || hour > 23) {
        std::string msg = "Invalid hour "+gammalib::str(hour)+" specified. "
                          "The hour needs to be comprised between 00 and 23";
        throw GException::invalid_argument(G_UTC, msg);
    }
    if (minute < 0 || minute > 59) {
        std::string msg = "Invalid minute "+gammalib::str(minute)+" specified. "
                          "The minute needs to be comprised between 00 and 59";
        throw GException::invalid_argument(G_UTC, msg);
    }
    if (second < 0 || second >= 60.0) {
        std::string msg = "Invalid second "+gammalib::str(second)+" specified. "
                          "The second needs to be comprised between 00 and <60";
        throw GException::invalid_argument(G_UTC, msg);
    }

    // Compute MJD day
    month--;
    for (int i = 0; i < month; ++i) { // Add days passed per month
        day += daymonth[i];
    }
    day += (year - 1972) * 365 - 1;   // Add days passed per year
    day += (year - 1969) / 4;         // Add leap days passed (every 4 years)
    day -= (year - 1901) / 100;       // Add leap days passed (not every 100 years)
    day += (year - 1601) / 400;       // Add leap days passed (every 400 years)
    day += 41317;                     // Add MJD at 1972

    // Compute MJD fraction (UTC)
    double fraction = ((double)hour * 3600.0 + (double)minute * 60.0 + second) *
                      gammalib::sec2day;

    // Compute MJD (UTC)
    double mjd = (double)day + fraction;

    // Get conversion from UTC to TAI to TT
    double correction = leap_seconds(mjd) + gammalib::tai2tt;

    // Convert MJD from UTC to TT
    mjd += correction * gammalib::sec2day;

    // Set time
    this->mjd(mjd);
    
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
 ***************************************************************************/
void GTime::set(const double& time, const GTimeReference& ref)
{
    // Convert time to specified time unit
    m_time = time * ref.unitseconds();
    
    // Compute time offset in seconds
    double offset = (mjd_ref - ref.mjdref()) * gammalib::sec_in_day;
        
    // Subtract time offset in seconds
    m_time -= offset;

    // Add leap seconds in case that time was given in UTC time system
    if (gammalib::toupper(ref.timesys()) == "UTC") {

        // Add leap seconds and offset between TAI and TT time system
        m_time += leap_seconds(this->mjd());
        m_time += gammalib::tai2tt;

    } // endif: Time was given in UTC time system

    // ... otherwise if system is TAI system then add TAI offset
    else if (gammalib::toupper(ref.timesys()) == "TAI") {
        m_time += gammalib::tai2tt;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set time from string
 *
 * @param[in] time Time string.
 * @param[in] ref Reference system.
 *
 * Sets the time from a string for a given reference system. The following
 * strings are valid:
 *
 *     "2016-10-05T15:08:56" (UTC string)
 *     "1800.0" (MET seconds in specified reference system)
 *     "1800.0 (TT)" (MET seconds in specified reference, TT system)
 *     "1800.0 (UTC)" (MET seconds in specified reference, UTC time system)
 *     "1800.0 (TAI)" (MET seconds in specified reference, TAI time system)
 *     "MJD 54609" (Modified Julian Days, TT system)
 *     "MJD 54609 (TT)" (Modified Julian Days, TT system)
 *     "MJD 54609 (UTC)" (Modified Julian Days, UTC system)
 *     "MJD 54609 (TAI)" (Modified Julian Days, TAI system)
 *     "JD 54609" (Julian Days, TT system)
 *     "JD 54609 (TT)" (Julian Days, TT system)
 *     "JD 54609 (UTC)" (Julian Days, UTC system)
 *     "JD 54609 (TAI)" (Julian Days, TAI system)
 *
 * If any other string is encountered, the numerical value is interpreted as
 * time is seconds using the TT system.
 *
 * Note that the TT, UTC or TAI attributes overwrite the values contained in
 * the specified reference system. The reference system is only used to
 * convert MET times in seconds.
 ***************************************************************************/
void GTime::set(const std::string& time, const GTimeReference& ref)
{
    // Strip any whitespace from string and convert it to upper case
    std::string str = gammalib::toupper(gammalib::strip_whitespace(time));

    // First check if the string is a UTC string
    long   year   = 0;
    int    month  = 0;
    long   day    = 0;
    long   hour   = 0;
    long   minute = 0;
    double second = 0.0;
    int    n      = sscanf(str.c_str(), "%ld-%d-%ldT%ld:%ld:%lg",
		                   &year, &month, &day, &hour, &minute, &second);
    if (n == 3 || n == 6) {
        utc(str);
    }

    // ... otherwise check for Modified Julian Days
    else if (str.find("MJD") == 0) {
        double      timeval = extract_timeval(str);
        std::string timesys = extract_timesys(str);
        if (timesys.empty()) {
            timesys = "TT";
        }
        mjd(timeval, timesys);
    }

    // ... otherwise check for Julian Days
    else if (str.find("JD") == 0) {
        double      timeval = extract_timeval(str);
        std::string timesys = extract_timesys(str);
        if (timesys.empty()) {
            timesys = "TT";
        }
        jd(timeval, timesys);
    }

    // ... otherwise take time as seconds and use the specified reference
    // system. If no TT, UTC or TAI arrtibutes are specified use the value
    // specified by the reference system.
    else {
        double      timeval = extract_timeval(str);
        std::string timesys = extract_timesys(str);
        if (timesys.empty()) {
            timesys = ref.timesys();
        }
        GTimeReference timeref(ref.mjdref(), ref.timeunit(), timesys, ref.timeref());
        set(timeval, timeref);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set time to current time
 *
 * Sets time to current time.
 ***************************************************************************/
void GTime::now(void)
{
    // Allocate variables
    struct std::tm timeStruct;
    std::time_t    now;
    char           buffer[100];

    // Get time
    now = std::time(NULL);
    #ifdef HAVE_GMTIME_R   
    std::gmtime_r(&now, &timeStruct);
    #else
    std::memcpy(&timeStruct, gmtime(&now), sizeof(struct tm));
    #endif

    // Write message type, time and task name to buffer
    std::sprintf(buffer, "%04d-%02d-%02dT%02d:%02d:%02d",
                         timeStruct.tm_year + 1900,
                         timeStruct.tm_mon + 1,
                         timeStruct.tm_mday,
                         timeStruct.tm_hour,
                         timeStruct.tm_min,
                         timeStruct.tm_sec);

    // Build string from buffer
    std::string date = buffer;

    // Set UTC time
    utc(date);

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
        result.append(gammalib::str(m_time)+" s (TT)");

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


/***********************************************************************//**
 * @brief Returns number of leap seconds for a given MJD
 *
 * @param[in] mjd Modified Julian Day in UTC time system.
 * @return Number of lead seconds.
 *
 * Return the number of leap seconds for a given MJD specified in the UTC
 * time system. This method returns valid number of leap seconds for the
 * years 1972-2017.
 *
 * See http://www.nist.gov/pml/div688/grp50/leapsecond.cfm for a table of
 * leap seconds.
 ***************************************************************************/
double GTime::leap_seconds(const double& mjd) const
{
    // Leap second table from 1972 on
    // see http://www.nist.gov/pml/div688/grp50/leapsecond.cfm
    // The first entry is 1-Jan-1972
    const long   leapsmjd[] = {41317,  // 1972-01-01
                               41498,  // 1972-06-30
                               41682,  // 1972-12-31
                               42047,  // 1973-12-31
                               42412,  // 1974-12-31
                               42777,  // 1975-12-31
                               43143,  // 1976-12-31
                               43508,  // 1977-12-31
                               43873,  // 1978-12-31
                               44238,  // 1979-12-31
                               44785,  // 1981-06-30
                               45150,  // 1982-06-30
                               45515,  // 1983-06-30
                               46246,  // 1985-06-30
                               47160,  // 1987-12-31
                               47891,  // 1989-12-31
                               48256,  // 1990-12-31
                               48803,  // 1992-06-30
                               49168,  // 1993-06-30
                               49533,  // 1994-06-30
                               50082,  // 1995-12-31
                               50629,  // 1997-06-30
                               51178,  // 1998-12-31
                               53735,  // 2005-12-31
                               54831,  // 2008-12-31
                               56108,  // 2012-06-30
                               57204,  // 2015-07-01
                               57754}; // 2017-01-01
    const double leapsecs[] = {10.0, 11.0, 12.0, 13.0, 14.0,
                               15.0, 16.0, 17.0, 18.0, 19.0,
                               20.0, 21.0, 22.0, 23.0, 24.0,
                               25.0, 26.0, 27.0, 28.0, 29.0,
                               30.0, 31.0, 32.0, 33.0, 34.0,
                               35.0, 36.0, 37.0};
    const int    n_leapsecs = sizeof(leapsmjd)/sizeof(long);

    // Extract MJD day and MJD fraction
    long day = (long)mjd;

    // Find the leap second MJD that is equal or larger to the specified
    // day
    int i = n_leapsecs - 1;
    while ((day <= leapsmjd[i]) && i > 0) {
        i--;
    }

    // Return leap seconds
    return (leapsecs[i]);
}


/***********************************************************************//**
 * @brief Extract time value from time string
 *
 * @param[in] time Time string.
 * @return Time value.
 *
 * Extracts the time value from a time string. The method strips any prefix
 * such as "MJD" or "JD" and any suffix starting with a left parentheses "("
 * and converts the remainder into a double precision value.
 ***************************************************************************/
double GTime::extract_timeval(const std::string& time) const
{
    // Find time
    size_t length = time.length();
    size_t start  = time.find_first_of("0123456789+-.");
    size_t stop   = time.find("(");

    // If there is a stop position then set length to stop
    if (stop != std::string::npos) {
        length = stop;
    }

    // Reduce length by start
    if (start != std::string::npos) {
        length -= start;
    }

    // Get substring containing the time value
    std::string value = time.substr(start, length);

    // Convert value into double precision
    double result = gammalib::todouble(value);

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Extract time system from time string
 *
 * @param[in] time Time string.
 * @return Time system.
 *
 * Extracts the time system from a time string. Valid time systems are:
 *
 *     "(TT)" (TT system)
 *     "(UTC)" (UTC system)
 *     "(TAI)" (TAI system)
 *
 * If no time system is found a blank string is returned.
 ***************************************************************************/
std::string GTime::extract_timesys(const std::string& time) const
{
    // Initialise time system with blank string
    std::string timesys = "";

    // Strip any whitespace from string and convert it to upper case
    std::string str = gammalib::toupper(gammalib::strip_whitespace(time));

    // Check for UTC
    if (str.find("(UTC)") != std::string::npos) {
        timesys = "UTC";
    }
    else if (str.find("(TAI)") != std::string::npos) {
        timesys = "TAI";
    }
    else if (str.find("(TT)") != std::string::npos) {
        timesys = "TT";
    }

    // Return time system
    return timesys;
}
