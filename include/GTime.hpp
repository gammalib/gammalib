/***************************************************************************
 *                          GTime.hpp - Time class                         *
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
    friend GTime operator+ (const GTime &a, const GTime &b);
    friend GTime operator- (const GTime &a, const GTime &b);
    friend GTime operator* (const double &a, const GTime &b);
    friend GTime operator* (const GTime &a, const double &b);
    friend GTime operator/ (const GTime &a, const double &b);
    friend bool  operator== (const GTime &a, const GTime &b);
    friend bool  operator!= (const GTime &a, const GTime &b);
    friend bool  operator< (const GTime &a, const GTime &b);
    friend bool  operator<= (const GTime &a, const GTime &b);
    friend bool  operator> (const GTime &a, const GTime &b);
    friend bool  operator>= (const GTime &a, const GTime &b);

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
    double         jd(void) const;
    double         mjd(void) const;
    double         secs(void) const;
    double         days(void) const;
    double         convert(const GTimeReference& ref) const;
    void           jd(const double& time);
    void           mjd(const double& time);
    void           secs(const double& seconds);
    void           days(const double& days);
    void           set(const double& time, const GTimeReference& ref);
    GTimeReference reference(void) const;
    std::string    print(void) const;
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GTime& time);
    void free_members(void);

    // Protected data members
    double m_time;   //!< Time in seconds in native reference
};


/***************************************************************************
 *                               Inline friends                            *
 ***************************************************************************/
inline
GTime operator+ (const GTime& a, const GTime& b)
{
    GTime result;
    result.m_time = a.m_time + b.m_time;
    return result;
}
inline
GTime operator- (const GTime& a, const GTime& b)
{
    GTime result;
    result.m_time = a.m_time - b.m_time;
    return result;
}
inline
GTime operator* (const double& a, const GTime& b)
{
    GTime result;
    result.m_time = a * b.m_time;
    return result;
}
inline
GTime operator* (const GTime& a, const double& b)
{
    GTime result;
    result.m_time = b * a.m_time;
    return result;
}
inline
GTime operator/ (const GTime& a, const double& b)
{
    GTime result;
    result.m_time = a.m_time / b;
    return result;
}
inline
bool operator== (const GTime &a, const GTime &b)
{
    return (a.m_time == b.m_time);
}
inline
bool operator!= (const GTime &a, const GTime &b)
{
    return (a.m_time != b.m_time);
}
inline
bool operator< (const GTime &a, const GTime &b)
{
    return (a.m_time < b.m_time);
}
inline
bool operator<= (const GTime &a, const GTime &b)
{
    return (a.m_time <= b.m_time);
}
inline
bool operator> (const GTime &a, const GTime &b)
{
    return (a.m_time > b.m_time);
}
inline
bool operator>= (const GTime &a, const GTime &b)
{
    return (a.m_time >= b.m_time);
}

#endif /* GTIME_HPP */
