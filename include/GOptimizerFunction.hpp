/***************************************************************************
 *    GOptimizerFunction.hpp  -  Optimizer function abstract base class    *
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
 * @file GOptimizerFunction.hpp
 * @brief Optimizer function abstract base class
 * @author J. Knodlseder
 */

#ifndef GOPTIMIZERFUNCTION_HPP
#define GOPTIMIZERFUNCTION_HPP

/* __ Includes ___________________________________________________________ */
#include "GOptimizerPars.hpp"
#include "GVector.hpp"
#include "GSparseMatrix.hpp"


/***********************************************************************//**
 * @class GOptimizerFunction
 *
 * @brief Optimizer function abstract base class
 *
 * This class provides an abstract interface for the function that is used
 * by the GOptimizer optimization class.
 * The method eval() returns the function value at a given set of parameters
 * that is defined by an instance of the optimizer parameter container class
 * GOptimizerPars.
 ***************************************************************************/
class GOptimizerFunction {

public:
    // Constructors and destructors
    GOptimizerFunction();
    GOptimizerFunction(const GOptimizerFunction& fct);
    virtual ~GOptimizerFunction();

    // Operators
    virtual GOptimizerFunction& operator= (const GOptimizerFunction& fct);

    // Virtual methods
    virtual void           eval(const GOptimizerPars& pars) = 0;
    virtual double*        value(void) = 0;
    virtual GVector*       gradient(void) = 0;
    virtual GSparseMatrix* covar(void) = 0;
 
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GOptimizerFunction& fct);
    void free_members(void);

};

#endif /* GOPTIMIZERFUNCTION_HPP */
