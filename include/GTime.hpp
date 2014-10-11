/***************************************************************************
 *                          GTime.hpp - Time class                         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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


/***********************************************************************//**
 * @class GTime
 *
 * @brief Handles times in a system independent way
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
    friend GTime operator+(const GTime &a, const GTime &b);
    friend GTime operator-(const GTime &a, const GTime &b);
    friend GTime operator*(const double &s, const GTime &time);
    friend GTime operator* (const GTime &time, const double &s);
    friend GTime operator/(const GTime &time, const double &s);
    friend bool  operator==(const GTime &a, const GTime &b);
    friend bool  operator!=(const GTime &a, const GTime &b);
    friend bool  operator<(const GTime &a, const GTime &b);
    friend bool  operator<=(const GTime &a, const GTime &b);
    friend bool  operator>(const GTime &a, const GTime &b);
    friend bool  operator>=(const GTime &a, const GTime &b);

public:
    // Constructors and destructors
    GTime(void);
    GTime(const GTime& time);
    explicit GTime(const double& time, const std::string& unit = "sec");
    explicit GTime(const double& time, const GTimeReference& ref);
    virtual ~GTime(void);
 
    // Operators
    GTime& operator=(const GTime& time);

    // Methods
    void           clear(void);
    GTime*         clone(void) const;
    std::string    classname(void) const;
    double         jd(void) const;
    double         mjd(void) const;
    const double&  secs(void) const;
    double         days(void) const;
    std::string    utc(void) const;
    double         convert(const GTimeReference& ref) const;
    void           jd(const double& time);
    void           mjd(const double& time);
    void           secs(const double& seconds);
    void           days(const double& days);
    void           utc(const std::string& time);
    void           set(const double& time, const GTimeReference& ref);
    GTimeReference reference(void) const;
    std::string    print(const GChatter& chatter = NORMAL) const;
  
protected:
    // Protected methods
    void   init_members(void);
    void   copy_members(const GTime& time);
    void   free_members(void);
    double leap_seconds(const double& mjd) const;

    // Protected data members
    double m_time;   //!< Time in seconds in native reference
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
 * @brief Return time in native reference (TT) (unit: seconds)
 *
 * @return Time in native reference [seconds].
 ***************************************************************************/
inline
const double& GTime::secs(void) const
{
    return m_time;
}


/***********************************************************************//**
 * @brief Set time in native reference in seconds (TT)
 *
 * @param[in] seconds Time (TT) [seconds].
 ***************************************************************************/
inline
void GTime::secs(const double& seconds)
{
    m_time = seconds;
    return;
}


/***********************************************************************//**
 * @brief Time addition operator friend
 *
 * @param[in] a First time.
 * @param[in] b Second time.
 * @return Added time.
 ***************************************************************************/
inline
GTime operator+(const GTime& a, const GTime& b)
{
    GTime result;
    result.m_time = a.m_time + b.m_time;
    return result;
}


/***********************************************************************//**
 * @brief Time subtraction operator friend
 *
 * @param[in] a First time.
 * @param[in] b Second time.
 * @return Subtracted time.
 ***************************************************************************/
inline
GTime operator-(const GTime& a, const GTime& b)
{
    GTime result;
    result.m_time = a.m_time - b.m_time;
    return result;
}


/***********************************************************************//**
 * @brief Time multiplication operator friend
 *
 * @param[in] s Multiplier.
 * @param[in] time Time.
 * @return Multiplied time.
 ***************************************************************************/
inline
GTime operator*(const double& s, const GTime& time)
{
    GTime result;
    result.m_time = s * time.m_time;
    return result;
}


/***********************************************************************//**
 * @brief Time multiplication operator friend
 *
 * @param[in] time Time.
 * @param[in] s Multiplier.
 * @return Multiplied time.
 ***************************************************************************/
inline
GTime operator* (const GTime& time, const double& s)
{
    GTime result;
    result.m_time = s * time.m_time;
    return result;
}


/***********************************************************************//**
 * @brief Time division operator friend
 *
 * @param[in] time Time.
 * @param[in] s Divider.
 * @return Divided time.
 ***************************************************************************/
inline
GTime operator/(const GTime& time, const double& s)
{
    GTime result;
    result.m_time = time.m_time / s;
    return result;
}


/***********************************************************************//**
 * @brief Time equality operator friend
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
 * @brief Time non-equality operator friend
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
 * @brief Time smaller than operator friend
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
 * @brief Time smaller than or equal operator friend
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
 * @brief Time larger than operator friend
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
 * @brief Time larger than or equal operator friend
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

#endif /* GTIME_HPP */
