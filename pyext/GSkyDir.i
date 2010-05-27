/***************************************************************************
 *               GSkyDir.i  -  Sky direction class SWIG file               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GSkyDir.i
 * @brief GSkyDir class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSkyDir.hpp"
%}


/***********************************************************************//**
 * @class GSkyDir
 *
 * @brief GSkyDir class interface defintion
 ***************************************************************************/
class GSkyDir {
public:
    // Constructors and destructors
    GSkyDir(void);
    GSkyDir(const GSkyDir& dir);
    virtual ~GSkyDir(void);

    // Methods
    void   clear(void);
    void   radec(const double& ra, const double& dec);
    void   radec_deg(const double& ra, const double& dec);
    void   lb(const double& l, const double& b);
    void   lb_deg(const double& l, const double& b);
    double l(void) const;
    double l_deg(void) const;
    double b(void) const;
    double b_deg(void) const;
    double ra(void) const;
    double ra_deg(void) const;
    double dec(void) const;
    double dec_deg(void) const;
    double dist(const GSkyDir& dir) const;
    double dist_deg(const GSkyDir& dir) const;
};


/***********************************************************************//**
 * @brief GSkyDir class extension
 ***************************************************************************/
%extend GSkyDir {
    char *__str__() {
        static char str_buffer[1001];
        std::ostringstream buffer;
        buffer << *self;
        std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 1001);
        str_buffer[1000] = '\0';
        return str_buffer;
    }
    GSkyDir copy() {
        return (*self);
    }
};
