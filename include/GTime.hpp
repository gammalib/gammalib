/***************************************************************************
 *                          GTime.hpp - Time class                         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2019 by Juergen Knoedlseder                         *
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
 * @file GTime.hpp
 * @brief Time class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GTIME_HPP
#define GTIME_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GTimeReference.hpp"

/* __ Forward declarations _______________________________________________ */


/***********************************************************************//**
 * @class GTime
 *
 * @brief Time class
 *
 * The GTime class stores a time value in seconds in a GammaLib native time
 * reference system. The GammaLib native time reference (i.e. time=0) is
 * defined as 
 *
 *                     January 1, 2010, 00:00:00 (TT)
 *
 * The time system is Terrestrial Time (TT). With respect to Coordinated
 * Universal Time (UTC), TT time is greater than UTC time by 66.184 sec at
 * January 1, 2010, 00:00:00. The difference is due to the introduction of
 * leap seconds that synchronize TT with the Earth rotation (UTC).
 ***************************************************************************/
class GTime : public GBase {

    // Operator friends
    friend GTime  operator+(const GTime&  time,    const double& seconds);
    friend GTime  operator+(const double& seconds, const GTime&  time);
    friend GTime  operator-(const GTime&  time,    const double& seconds);
    friend double operator-(const GTime&  a,       const GTime&  b);
    friend bool   operator==(const GTime& a,       const GTime& b);
    friend bool   operator!=(const GTime& a,       const GTime& b);
    friend bool   operator<(const GTime&  a,       const GTime& b);
    friend bool   operator<=(const GTime& a,       const GTime& b);
    friend bool   operator>(const GTime&  a,       const GTime& b);
    friend bool   operator>=(const GTime& a,       const GTime& b);

public:
    // Constructors and destructors
    GTime(void);
    GTime(const GTime& time);
    GTime(const double& time, const std::string& unit = "sec");
    GTime(const double& time, const GTimeReference& ref);
    GTime(const std::string& time, const GTimeReference& ref);
    explicit GTime(const std::string& time);
    virtual ~GTime(void);

    // Operators
    GTime& operator=(const GTime& time);
    GTime& operator+=(const double& seconds);
    GTime& operator-=(const double& seconds);

    // Methods
    void           clear(void);
    GTime*         clone(void) const;
    std::string    classname(void) const;
    double         jd(void) const;
    double         jd(const std::string& timesys) const;
    double         mjd(void) const;
    double         mjd(const std::string& timesys) const;
    const double&  secs(void) const;
    double         secs(const std::string& timesys) const;
    double         days(void) const;
    double         days(const std::string& timesys) const;
    std::string    utc(const int& precision = 0) const;
    double         gmst(void) const;
    double         gast(void) const;
    double         lmst(const double& geolon) const;
    double         last(const double& geolon) const;
    double         convert(const GTimeReference& ref) const;
    void           jd(const double& time);
    void           jd(const double& time, const std::string& timesys);
    void           mjd(const double& time);
    void           mjd(const double& time, const std::string& timesys);
    void           secs(const double& seconds);
    void           secs(const double& seconds, const std::string& timesys);
    void           days(const double& days);
    void           days(const double& days, const std::string& timesys);
    void           utc(const std::string& time);
    void           set(const double& time, const GTimeReference& ref);
    void           set(const std::string& time, const GTimeReference& ref);
    void           set(const std::string& time);
    void           now(void);
    GTimeReference reference(void) const;
    std::string    print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const GTime& time);
    void        free_members(void);
    double      leap_seconds(const double& mjd) const;
    bool        is_leap_year(const int& year) const;
    int         days_in_year(const int& year) const;
    double      extract_timeval(const std::string& time) const;
    std::string extract_timesys(const std::string& time) const;

    // Protected data members
    double m_time; //!< Time in seconds in native reference (TT)
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GTime").
 ***************************************************************************/
inline
std::string GTime::classname(void) const
{
    return ("GTime");
}


/***********************************************************************//**
 * @brief Return time in seconds in native reference (TT)
 *
 * @return Time in native reference (seconds).
 ***************************************************************************/
inline
const double& GTime::secs(void) const
{
    return m_time;
}


/***********************************************************************//**
 * @brief Set time in seconds in native reference (TT)
 *
 * @param[in] seconds Time in native reference (seconds).
 ***************************************************************************/
inline
void GTime::secs(const double& seconds)
{
    m_time = seconds;
    return;
}


/***********************************************************************//**
 * @brief Signals if year is a leap year
 *
 * @param[in] year Year (four digits integer).
 * @return True if year is a leap year.
 ***************************************************************************/
inline
bool GTime::is_leap_year(const int& year) const
{
    bool is_leap_year = (((year % 400) == 0) || 
                         (((year % 100) != 0) && ((year % 4) == 0)));
    return (is_leap_year);
}


/***********************************************************************//**
 * @brief Returns number of days in year
 *
 * @param[in] year Year (four digits integer).
 * @return Number of days in year.
 ***************************************************************************/
inline
int GTime::days_in_year(const int& year) const
{
    int days = (is_leap_year(year)) ? 366 : 365;
    return (days);
}


/***********************************************************************//**
 * @brief Add seconds to time
 *
 * @param[in] seconds Seconds.
 * @return Time.
 *
 * Adds @p seconds to the time.
 ***************************************************************************/
inline
GTime& GTime::operator+=(const double& seconds)
{
    m_time += seconds;
    return *this;
}


/***********************************************************************//**
 * @brief Subtract seconds from time
 *
 * @param[in] seconds Seconds.
 * @return Time.
 *
 * Subtracts @p seconds from the time.
 ***************************************************************************/
inline
GTime& GTime::operator-=(const double& seconds)
{
    m_time -= seconds;
    return *this;
}


/***********************************************************************//**
 * @brief Add seconds to time
 *
 * @param[in] time Time.
 * @param[in] seconds Seconds.
 * @return Time.
 *
 * Adds @p seconds to the @p time.
 ***************************************************************************/
inline
GTime operator+(const GTime& time, const double& seconds)
{
    GTime result;
    result.m_time = time.m_time + seconds;
    return result;
}


/***********************************************************************//**
 * @brief Add seconds to time
 *
 * @param[in] seconds Seconds.
 * @param[in] time Time.
 * @return Time.
 *
 * Adds @p seconds to the @p time.
 ***************************************************************************/
inline
GTime operator+(const double& seconds, const GTime& time)
{
    GTime result;
    result.m_time = time.m_time + seconds;
    return result;
}


/***********************************************************************//**
 * @brief Subtract seconds from time
 *
 * @param[in] time Time.
 * @param[in] seconds Seconds.
 * @return Time.
 *
 * Subtracts @p seconds from the @p time.
 ***************************************************************************/
inline
GTime operator-(const GTime& time, const double& seconds)
{
    GTime result;
    result.m_time = time.m_time - seconds;
    return result;
}


/***********************************************************************//**
 * @brief Subtract times
 *
 * @param[in] a First time.
 * @param[in] b Second time.
 * @return Time difference in seconds.
 *
 * Subtracts time @p b from time @p a and returns the difference in seconds.
 ***************************************************************************/
inline
double operator-(const GTime& a, const GTime& b)
{
    return (a.m_time - b.m_time);
}
 

/***********************************************************************//**
 * @brief Check if times are equal
 *
 * @param[in] a First time.
 * @param[in] b Second time.
 * @return True if @p a equals @p b.
 ***************************************************************************/
inline
bool operator==(const GTime &a, const GTime &b)
{
    return (a.m_time == b.m_time);
}


/***********************************************************************//**
 * @brief Check if times are not equal
 *
 * @param[in] a First time.
 * @param[in] b Second time.
 * @return True if @p a is not equal to @p b.
 ***************************************************************************/
inline
bool operator!=(const GTime &a, const GTime &b)
{
    return (a.m_time != b.m_time);
}


/***********************************************************************//**
 * @brief Check if time is smaller than other time
 *
 * @param[in] a First time.
 * @param[in] b Second time.
 * @return True if @p a is smaller than @p b.
 ***************************************************************************/
inline
bool operator<(const GTime &a, const GTime &b)
{
    return (a.m_time < b.m_time);
}


/***********************************************************************//**
 * @brief Check if time is smaller than or equal to other time
 *
 * @param[in] a First time.
 * @param[in] b Second time.
 * @return True if @p a is smaller than or equal to @p b.
 ***************************************************************************/
inline
bool operator<=(const GTime &a, const GTime &b)
{
    return (a.m_time <= b.m_time);
}


/***********************************************************************//**
 * @brief Check if time is larger than other time
 *
 * @param[in] a First time.
 * @param[in] b Second time.
 * @return True if @p a is larger than @p b.
 ***************************************************************************/
inline
bool operator>(const GTime &a, const GTime &b)
{
    return (a.m_time > b.m_time);
}


/***********************************************************************//**
 * @brief Check if time is larger than or equal to other time
 *
 * @param[in] a First time.
 * @param[in] b Second time.
 * @return True if @p a is larger than or equal to @p b.
 ***************************************************************************/
inline
bool operator>=(const GTime &a, const GTime &b)
{
    return (a.m_time >= b.m_time);
}


/***********************************************************************//**
 * @brief Set time for string for native time reference system
 *
 * @param[in] time Time string.
 ***************************************************************************/
inline
void GTime::set(const std::string& time)
{
    set(time, this->reference());
    return;
}

#endif /* GTIME_HPP */
