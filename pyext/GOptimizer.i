/***************************************************************************
 *            GOptimizer.i  -  Optimizer class SWIG interface              *
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
 * @file GOptimizer.i
 * @brief GOptimizer class SWIG interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GOptimizer.hpp"
%}


/***********************************************************************//**
 * @class GOptimizer
 *
 * @brief GOptimizer class SWIG interface defintion.
 ***************************************************************************/
class GOptimizer {
public:
    // Constructors and destructors
    GOptimizer();
    GOptimizer(const GOptimizer& opt);
    virtual ~GOptimizer();

    // Operators (we need those to make the class abstract!!!)
    virtual GOptimizerPars& operator() (GOptimizerFunction& fct, GOptimizerPars& p) = 0;
    virtual GModels&        operator() (GOptimizerFunction& fct, GModels& m) = 0;
};
