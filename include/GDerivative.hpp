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
#include <string>
#include <iostream>
#include "GLog.hpp"
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

    // I/O friends
    friend std::ostream& operator<<(std::ostream& os, const GDerivative& dx);
    friend GLog&         operator<<(GLog& os,         const GDerivative& dx);

public:

    // Constructors and destructors
    GDerivative(void);
    explicit GDerivative(GFunction* func);
    GDerivative(const GDerivative& dx);
    virtual ~GDerivative(void);

    // Operators
    GDerivative& operator= (const GDerivative& dx);

    // Methods
    void             clear(void);
    GDerivative*     clone(void) const;
    void             max_iter(const int& max_iter) { m_max_iter=max_iter; }
    void             eps(const double& eps) { m_eps=eps; }
    void             step_frac(const double& f) { m_step_frac=f; }
    void             silent(const bool& silent) { m_silent=silent; }
    const int&       iter(void) const { return m_iter; }
    const int&       max_iter(void) const { return m_max_iter; }
    const double&    eps(void) const { return m_eps; }
    const double&    step_frac(void) const { return m_step_frac; }
    const bool&      silent(void) const { return m_silent; }
    void             function(GFunction* func) { m_func = func; }
    const GFunction* function(void) const { return m_func; }
    double           value(const double& x, double step = 0.0);
    double           ridder(const double& x, const double& h, double& err);
    double           minuit(const double& x, double& err);
    std::string      print(void) const;

protected:
    // Protected methods
    void   init_members(void);
    void   copy_members(const GDerivative& dx);
    void   free_members(void);

    // Protected members
    GFunction* m_func;         //!< Pointer to function
    double     m_eps;          //!< Derivative precision
    double     m_step_frac;    //!< Value fraction to use for initial step
    int        m_max_iter;     //!< Maximum number of iterations
    int        m_iter;         //!< Number of iterations used
    bool       m_silent;       //!< Suppress warnings
};

#endif /* GDERIVATIVE_HPP */
