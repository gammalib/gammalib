/***************************************************************************
 *      GOptimizerPars.i  -  Parameter container class SWIG interface      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GOptimizerPars.i
 * @brief GOptimizerPars class SWIG interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GOptimizerPars.hpp"
%}


/***********************************************************************//**
 * @class GOptimizerPars
 *
 * @brief GOptimizerPars container class SWIG interface defintion.
 ***************************************************************************/
class GOptimizerPars {
public:
    // Constructors and destructors
    GOptimizerPars(void);
    GOptimizerPars(const GOptimizerPars& pars);
    ~GOptimizerPars();

    // Methods
    int        npars(void) const { return m_npars; }
    int        nfree(void) const;
    GModelPar* par(int index) const;
};


/***********************************************************************//**
 * @brief GOptimizerPars class extension
 ***************************************************************************/
%extend GOptimizerPars {
    char *__str__() {
        static char str_buffer[10001];
        std::ostringstream buffer;
        buffer << *self;
        std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 10001);
        str_buffer[10000] = '\0';
        return str_buffer;
    }
    GOptimizerPars copy() {
        return (*self);
    }
};
