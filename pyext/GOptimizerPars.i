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
 * @brief Optimizer parameter container class Python interface
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
 * @brief Optimizer parameter container class
 ***************************************************************************/
class GOptimizerPars {
public:
    // Constructors and destructors
    GOptimizerPars(void);
    GOptimizerPars(const GOptimizerPars& pars);
    virtual ~GOptimizerPars(void);

    // Methods
    int              npars(void) const;
    int              nfree(void) const;
    GModelPar&       par(const int& index);
    //const GModelPar& par(const int& index) const; // now overloaded by SWIG
};


/***********************************************************************//**
 * @brief GOptimizerPars class extension
 ***************************************************************************/
%extend GOptimizerPars {
    char *__str__() {
        return tochar(self->print());
    }
    GOptimizerPars copy() {
        return (*self);
    }
};
