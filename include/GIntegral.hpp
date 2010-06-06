/***************************************************************************
 *                  GIntegral.hpp  -  Integration class                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
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
 * @brief GIntegral class interface definition.
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
    void              integrand(GIntegrand* integrand) { m_integrand = integrand; }
    const GIntegrand* integrand(void) const { return m_integrand; }
    double            romb(double a, double b);

protected:
    // Protected methods
    void   init_members(void);
    void   copy_members(const GIntegral& integral);
    void   free_members(void);
    void   polint(double* xa, double* ya, int n, double x, double *y, double *dy);
    double trapzd(double a, double b, int n);
    
    // Protected data area
    GIntegrand* m_integrand;    //!< Pointer to integrand
    double      m_eps;          //!< Integration precision
    int         m_max_iter;     //!< Maximum number of iterations

};

#endif /* GINTEGRAL_HPP */
