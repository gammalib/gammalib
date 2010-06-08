/***************************************************************************
 *                          GTime.hpp - Time class                         *
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
 * @file GTime.hpp
 * @brief Time value class definition.
 * @author J. Knodlseder
 */

#ifndef GTIME_HPP
#define GTIME_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>


/***********************************************************************//**
 * @class GTime
 *
 * @brief Class that handles times in a system independent way.
 *
 * The GTime class stores a time value in MJD and implements methods that
 * provide automatic conversion of the time values in other systems. This
 * makes instrument specific implementations more robust and reduces the
 * risk of unit errors.
 ***************************************************************************/
class GTime {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GTime& time);

    // Operator friends
    friend GTime operator+ (const GTime &a, const GTime &b);
    friend GTime operator- (const GTime &a, const GTime &b);
    friend GTime operator* (const double &a, const GTime &b);
    friend GTime operator* (const GTime &a, const double &b);
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
    void   clear(void) { m_time = 0.0; }
    double mjd(void) const;
    void   mjd(const double& time);
    double met(void) const;
    void   met(const double& time);
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GTime& time);
    void free_members(void);

    // Protected data members
    double m_time;          //!< Time in MJD
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
