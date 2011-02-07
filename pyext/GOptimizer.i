/***************************************************************************
 *            GOptimizer.i  -  Optimizer class Python interface            *
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
 * @file GOptimizer.i
 * @brief GOptimizer class Python interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GOptimizer.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GOptimizer
 *
 * @brief GOptimizer class Python interface defintion
 ***************************************************************************/
class GOptimizer {
public:
    // Constructors and destructors
    GOptimizer(void);
    GOptimizer(const GOptimizer& opt);
    virtual ~GOptimizer();

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GOptimizer* clone(void) const = 0;
    virtual std::string print(void) const = 0;
};


/***********************************************************************//**
 * @brief GOptimizer class extension
 ***************************************************************************/
%extend GOptimizer {
    char *__str__() {
        return tochar(self->print());
    }
};
