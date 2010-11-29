/***************************************************************************
 *                     GTime.i  -  Time class python I/F                   *
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
 * @file GTime.i
 * @brief GTime class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GTime.hpp"
%}


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
    // Operator friends
    /*
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
    */

public:
    // Constructors and destructors
    GTime(void);
    GTime(const GTime& time);
    ~GTime(void);
 
    // Methods
    void   clear(void) { m_time = 0.0; }
    double mjd(void) const;
    void   mjd(const double& time);
    double met(void) const;
    void   met(const double& time);
};


/***********************************************************************//**
 * @brief GTime class extension
 ***************************************************************************/
%extend GTime {
    /*
    char *__str__() {
        static std::string result = self->print();
        return ((char*)result.c_str());
    }
    */
    GTime copy() {
        return (*self);
    }
};
