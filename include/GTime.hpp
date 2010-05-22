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

public:
    // Constructors and destructors
    GTime(void);
    GTime(const GTime& time);
    ~GTime(void);
 
    // Operators
    GTime& operator= (const GTime& time);

    // Methods
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

#endif /* GTIME_HPP */
