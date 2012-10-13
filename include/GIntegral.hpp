/***************************************************************************
 *                  GIntegral.hpp  -  Integration class                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
 * @file GIntegral.hpp
 * @brief Integration class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GINTEGRAL_HPP
#define GINTEGRAL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GBase.hpp"
#include "GLog.hpp"
#include "GIntegrand.hpp"


/***********************************************************************//**
 * @class GIntegral
 *
 * @brief GIntegral class interface defintion.
 *
 * This class allows to perform integration using various methods. The
 * integrand is implemented by a derived class of GIntegrand.
 ***************************************************************************/
class GIntegral : public GBase {

    // I/O friends
    friend std::ostream& operator<<(std::ostream& os, const GIntegral& integral);
    friend GLog&         operator<<(GLog& log,        const GIntegral& integral);

public:

    // Constructors and destructors
    explicit GIntegral(void);
    explicit GIntegral(GIntegrand* integrand);
    GIntegral(const GIntegral& integral);
    virtual ~GIntegral(void);

    // Operators
    GIntegral& operator=(const GIntegral& integral);

    // Methods
    void              clear(void);
    GIntegral*        clone(void) const;
    void              max_iter(const int& max_iter) { m_max_iter=max_iter; }
    void              eps(const double& eps) { m_eps=eps; }
    void              silent(const bool& silent) { m_silent=silent; }
    const int&        iter(void) const { return m_iter; }
    const int&        max_iter(void) const { return m_max_iter; }
    const double&     eps(void) const { return m_eps; }
    const bool&       silent(void) const { return m_silent; }
    void              integrand(GIntegrand* integrand) { m_integrand=integrand; }
    const GIntegrand* integrand(void) const { return m_integrand; }
    double            romb(double a, double b, int k = 5);
    double            trapzd(double a, double b, int n = 1, double result = 0.0);
    std::string       print(void) const;

protected:
    // Protected methods
    void   init_members(void);
    void   copy_members(const GIntegral& integral);
    void   free_members(void);
    double polint(double* xa, double* ya, int n, double x, double *dy);

    // Protected data area
    GIntegrand* m_integrand;    //!< Pointer to integrand
    double      m_eps;          //!< Integration precision
    int         m_max_iter;     //!< Maximum number of iterations
    int         m_iter;         //!< Number of iterations used
    bool        m_silent;       //!< Suppress integration warnings
};

#endif /* GINTEGRAL_HPP */
