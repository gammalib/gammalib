/***************************************************************************
 *           GRan.i  -  Random number generator class python I/F           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GRan.i
 * @brief GRan class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GRan.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GRan
 *
 * @brief Random number generator class
 *
 * The GRan class implements a random number generator that is inspired from
 * the Ran structure given in Numerical Recipes, Third Edition (p. 341ff).
 ***************************************************************************/
class GRan {
public:
    // Constructors and destructors
    GRan(void);
    GRan(unsigned long long int seed);
    GRan(const GRan& ran);
    virtual ~GRan(void);
 
    // Methods
    void                   clear(void);
    GRan*                  clone(void) const;
    void                   seed(unsigned long long int seed);
    unsigned long int      int32(void);
    unsigned long long int int64(void);
    double                 uniform(void);
    double                 exp(const double& lambda);
};


/***********************************************************************//**
 * @brief GRan class extension
 ***************************************************************************/
%extend GRan {
    char *__str__() {
        return tochar(self->print());
    }
    GRan copy() {
        return (*self);
    }
};
