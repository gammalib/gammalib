/***************************************************************************
 *                  GDerivative.hpp  -  Derivative class                   *
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
 * @file GDerivative.hpp
 * @brief GDerivative class interface definition.
 * @author J. Knodlseder
 */

#ifndef GDERIVATIVE_HPP
#define GDERIVATIVE_HPP

/* __ Includes ___________________________________________________________ */
#include "GFunction.hpp"


/***********************************************************************//**
 * @class GDerivative
 *
 * @brief Numerical derivatives class interface defintion
 *
 * This class allows to compute numerical derivatives using various methods.
 * The function to be derived is implemented by the abstract GFunction
 * class.
 ***************************************************************************/
class GDerivative {

public:

    // Constructors and destructors
    GDerivative(void);
    explicit GDerivative(GFunction* func);
    GDerivative(const GDerivative& dx);
    virtual ~GDerivative(void);

    // Operators
    GDerivative& operator= (const GDerivative& dx);

    // Methods
    void             function(GFunction* func) { m_func = func; }
    const GFunction* function(void) const { return m_func; }
    double           value(const double& x);
    double           ridder(const double& x, const double& h, double& err);

protected:
    // Protected methods
    void   init_members(void);
    void   copy_members(const GDerivative& dx);
    void   free_members(void);

    // Protected data area
    GFunction* m_func;    //!< Pointer to function
};

#endif /* GDERIVATIVE_HPP */
