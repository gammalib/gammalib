/***************************************************************************
 *                  GIntegral.hpp  -  Integration class                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GIntegral.hpp
 * @brief Integration class interface definition
 * @author J. Knodlseder
 */

#ifndef GINTEGRAL_HPP
#define GINTEGRAL_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GIntegrand.hpp"


/***********************************************************************//**
 * @class GIntegral
 *
 * @brief GIntegral class interface defintion.
 *
 * This class allows to perform integration using various methods. The
 * integrand is implemented by a derived class of GIntegrand.
 ***************************************************************************/
class GIntegral {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GIntegral& integral);

public:

    // Constructors and destructors
    explicit GIntegral(void);
    explicit GIntegral(GIntegrand* integrand);
    GIntegral(const GIntegral& integral);
    virtual ~GIntegral(void);

    // Operators
    GIntegral& operator= (const GIntegral& integral);

    // Methods
    void              max_iter(const int& max_iter) { m_max_iter=max_iter; }
    void              eps(const double& eps) { m_eps=eps; }
    int               max_iter(void) const { return m_max_iter; }
    double            eps(void) const { return m_eps; }
    int               iter(void) const { return m_iter; }
    void              integrand(GIntegrand* integrand) { m_integrand = integrand; }
    const GIntegrand* integrand(void) const { return m_integrand; }
    double            romb(double a, double b, int k = 5);

protected:
    // Protected methods
    void   init_members(void);
    void   copy_members(const GIntegral& integral);
    void   free_members(void);
    double polint(double* xa, double* ya, int n, double x, double *dy);
    double trapzd(double a, double b, int n, double result);

    // Protected data area
    GIntegrand* m_integrand;    //!< Pointer to integrand
    double      m_eps;          //!< Integration precision
    int         m_max_iter;     //!< Maximum number of iterations
    int         m_iter;         //!< Number of iterations used
};

#endif /* GINTEGRAL_HPP */
