/***************************************************************************
 *  GOptimizerFunction.hpp  -  Abstract base class for optimizer function  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009 by Jurgen Knodlseder                                *
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
 * @brief GOptimizerFunction abstract base class interface definition.
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
 * @brief GOptimizerFunction abstract base class interface defintion.
 ***************************************************************************/
class GOptimizerFunction {

public:
    // Constructors and destructors
    GOptimizerFunction(const GOptimizerPars& pars);
    GOptimizerFunction(const GOptimizerFunction& fct);
    virtual ~GOptimizerFunction();

    // Operators
    virtual GOptimizerFunction& operator= (const GOptimizerFunction& fct);

    // Virtual methods
    virtual double        value(void) const = 0;
    virtual GVector       gradient(void) const = 0;
    virtual GSparseMatrix covar(void) const = 0;
 
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GOptimizerFunction& fct);
    void free_members(void);

};

#endif /* GOPTIMIZERFUNCTION_HPP */
