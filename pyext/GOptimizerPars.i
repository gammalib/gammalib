/***************************************************************************
 *      GOptimizerPars.i  -  Parameter container class SWIG interface      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2011 by Jurgen Knodlseder                           *
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
#include "GTools.hpp"
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
        std::ostringstream buffer;
        buffer << *self;
        std::string str = buffer.str();
        return tochar(str);
    }
    GOptimizerPars copy() {
        return (*self);
    }
};
