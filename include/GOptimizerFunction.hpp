/***************************************************************************
 *     GOptimizerFunction.hpp - Optimizer function abstract base class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2013 by Juergen Knoedlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GOptimizerFunction.hpp
 * @brief Optimizer function abstract base class
 * @author Juergen Knoedlseder
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
 *
 * The eval() method returns the function value at a given set of parameters
 * that is defined by an instance of the optimizer parameter container class
 * GOptimizerPars. The value() method returns the actual function value at
 * these parameters, and the gradient() and covar() methods return pointers
 * on the gradient vector and the covariance matrix at the parameter values.
 ***************************************************************************/
class GOptimizerFunction {

public:
    // Constructors and destructors
    GOptimizerFunction(void);
    GOptimizerFunction(const GOptimizerFunction& fct);
    virtual ~GOptimizerFunction(void);

    // Operators
    virtual GOptimizerFunction& operator= (const GOptimizerFunction& fct);

    // Virtual methods
    virtual void           eval(const GOptimizerPars& pars) = 0;
    virtual double         value(void) = 0;
    virtual GVector*       gradient(void) = 0;
    virtual GSparseMatrix* covar(void) = 0;
 
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GOptimizerFunction& fct);
    void free_members(void);
};

#endif /* GOPTIMIZERFUNCTION_HPP */
