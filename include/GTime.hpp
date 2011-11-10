/***************************************************************************
 *                          GTime.hpp - Time class                         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Juergen Knoedlseder                         *
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
 * @author J. Knoedlseder
 */

#ifndef GTIME_HPP
#define GTIME_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include <iostream>
#include "GLog.hpp"


/***********************************************************************//**
 * @class GTime
 *
 * @brief Class that handles times in a system independent way
 *
 * The GTime class stores a time value in seconds and its MJD reference in
 * days. The time can be retrieved in any MJD reference.
 *
 * @todo Make use of MJD reference in operators.
 ***************************************************************************/
class GTime {

    // I/O friends
    friend GLog&         operator<< (GLog& log,        const GTime& time);
    friend std::ostream& operator<< (std::ostream& os, const GTime& time);

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
    explicit GTime(const double&      time,
                   const double&      mrdref,
                   const std::string& timeunit,
                   const std::string& timesys = "TT",
                   const std::string& timeref = "local");
    explicit GTime(const double&      time,
                   const int&         mjdrefi,
                   const double&      mrdreff,
                   const std::string& timeunit,
                   const std::string& timesys = "TT",
                   const std::string& timeref = "local");
    virtual ~GTime(void);
 
    // Operators
    GTime& operator= (const GTime& time);

    // Methods
    void        clear(void);
    double      jd(void) const;
    double      mjd(void) const;
    double      met(void) const;
    void        jd(const double& time);
    void        mjd(const double& time);
    void        met(const double& time);
    void        time(const double& time,
                     const double&      mrdref,
                     const std::string& timeunit,
                     const std::string& timesys = "TT",
                     const std::string& timeref = "local");
    void        time(const double&      time,
                     const int&         mjdrefi,
                     const double&      mrdreff,
                     const std::string& timeunit,
                     const std::string& timesys = "TT",
                     const std::string& timeref = "local");
    double      time(void) const;
    double      time(const double& mrdref) const;
    double      time(const int& mjdrefi, const double& mrdreff) const;
    double      mjdref(void) const;
    std::string print(void) const;
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GTime& time);
    void free_members(void);

    // Protected data members
    double m_time;          //!< Time in seconds
    double m_mjdref;        //!< Time reference
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


/***************************************************************************
 *                                 Typedefs                                *
 ***************************************************************************/
typedef std::vector<GTime> GTimes;

#endif /* GTIME_HPP */
