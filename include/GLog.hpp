/***************************************************************************
 *                       GLog.hpp - Information logger                     *
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
 * @file GLog.hpp
 * @brief Information logger class definition
 * @author Jurgen Knodlseder
 */

#ifndef GLOG_HPP
#define GLOG_HPP

/* __ Includes ___________________________________________________________ */
#include <time.h>
#include <vector>
#include <string>
#include <iostream>


/***********************************************************************//**
 * @class GLog
 *
 * @brief Information logger interface defintion.
 ***************************************************************************/
class GLog {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GLog& log);

public:
    // Constructors and destructors
    GLog(void);
    GLog(const GLog& log);
    ~GLog(void);
 
    // Operators
    GLog& operator= (const GLog& log);

    // Methods
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GLog& log);
    void free_members(void);

    // Protected data members
};

#endif /* GLOG_HPP */
