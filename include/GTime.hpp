/***************************************************************************
 *                          GTime.hpp - Time class                         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GTime.hpp
 * @brief Time value class definition.
 * @author J. Knodlseder
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
 * @brief Class that handles times in a system independent way.
 *
 * The GTime class stores a time value in MET which are defined as the
 * number of seconds since 2001-01-01 00:00:00.000 UTC in TT (Terrestrial
 * Time). This is the Fermi/LAT reference time. See
 * http://fermi.gsfc.nasa.gov/ssc/data/analysis/documentation/Cicerone/Cicerone_Data/Time_in_ScienceTools.html
 * GTime implements methods that provide automatic conversion of the
 * time values in other systems. This makes instrument specific
 * implementations more robust and reduces the risk of unit errors.
 *
 * @todo Add general conversion routine that makes use of standard
 *       HEASARC keyword values, such as TIMESYS, MJDREF (or MJDREFI/MJDREFF),
 *       TIMEUNIT, TIMEZERO, TIMEREF.
 ***************************************************************************/
class GTime {

    // I/O friends
    friend GLog&         operator<< (GLog& log, const GTime& time);
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
    ~GTime(void);
 
    // Operators
    GTime& operator= (const GTime& time);

    // Methods
    void        clear(void) { m_time = 0.0; }
    double      jd(void) const;
    double      mjd(void) const;
    double      met(void) const;
    void        jd(const double& time);
    void        mjd(const double& time);
    void        met(const double& time);
    std::string print(void) const;
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GTime& time);
    void free_members(void);

    // Protected data members
    double m_time;          //!< Time in MET
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
