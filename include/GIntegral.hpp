/***************************************************************************
 *                   GIntegral.hpp - Integration class                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
#include <vector>
#include "GBase.hpp"
#include "GFunction.hpp"


/***********************************************************************//**
 * @class GIntegral
 *
 * @brief GIntegral class interface definition.
 *
 * This class allows to perform integration using various methods. The
 * integrand is implemented by a derived class of GFunction.
 ***************************************************************************/
class GIntegral : public GBase {

public:

    // Constructors and destructors
    explicit GIntegral(void);
    explicit GIntegral(GFunction* kernel);
    GIntegral(const GIntegral& integral);
    virtual ~GIntegral(void);

    // Operators
    GIntegral& operator=(const GIntegral& integral);

    // Methods
    void               clear(void);
    GIntegral*         clone(void) const;
    std::string        classname(void) const;
    void               max_iter(const int& iter);
    const int&         max_iter(void) const;
    void               fixed_iter(const int& iter);
    const int&         fixed_iter(void) const;
    void               eps(const double& eps);
    const double&      eps(void) const;
    void               silent(const bool& silent);
    const bool&        silent(void) const;
    const int&         iter(void) const;
    const int&         calls(void) const;
    const bool&        is_valid(void) const;
    const std::string& message(void) const;
    void               kernel(GFunction* kernel);
    const GFunction*   kernel(void) const;
    double             romberg(std::vector<double> bounds,
                               const int& order = 5);
    double             romberg(const double& a, const double& b,
                               const int& order = 5);
    double             trapzd(const double& a, const double& b,
                              const int& n = 1, double result = 0.0);
    double             adaptive_simpson(const double& a, const double& b) const;
    double             gauss_kronrod(const double& a, const double& b) const;
    std::string        print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void   init_members(void);
    void   copy_members(const GIntegral& integral);
    void   free_members(void);
    double polint(double* xa, double* ya, int n, double x, double *dy);
    double adaptive_simpson_aux(const double& a, const double& b,
                                const double& eps, const double& S,
                                const double& fa, const double& fb,
                                const double& fc,
                                const int& bottom) const;
    double rescale_error(double err,
                         const double& result_abs,
                         const double& result_asc) const;

    // Protected data area
    GFunction*  m_kernel;    //!< Pointer to function kernel
    double      m_eps;       //!< Requested relative integration precision
    int         m_max_iter;  //!< Maximum number of iterations
    int         m_fix_iter;  //!< Fixed number of iterations
    bool        m_silent;    //!< Suppress integration warnings in console

    // Integrator results
    mutable int         m_iter;       //!< Number of iterations used
    mutable int         m_calls;      //!< Number of function calls used
    mutable bool        m_isvalid;    //!< Integration result valid (true=yes)
    mutable bool        m_has_abserr; //!< Has absolute integration error
    mutable bool        m_has_relerr; //!< Has relative integration error
    mutable double      m_abserr;     //!< Absolute integration error
    mutable double      m_relerr;     //!< Absolute integration error
    mutable std::string m_message;    //!< Status message (if result is invalid)

};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GIntegral").
 ***************************************************************************/
inline
std::string GIntegral::classname(void) const
{
    return ("GIntegral");
}


/***********************************************************************//**
 * @brief Return number of iterations
 *
 * @return Number of iterations.
 ***************************************************************************/
inline
const int& GIntegral::iter(void) const
{
    return m_iter;
}


/***********************************************************************//**
 * @brief Set maximum number of iterations
 *
 * @param[in] iter Maximum number of iterations.
 ***************************************************************************/
inline
void GIntegral::max_iter(const int& iter)
{
    m_max_iter = iter;
    return;
}


/***********************************************************************//**
 * @brief Return maximum number of iterations
 *
 * @return Maximum number of iterations.
 ***************************************************************************/
inline
const int& GIntegral::max_iter(void) const
{
    return m_max_iter;
}


/***********************************************************************//**
 * @brief Set fixed number of iterations
 *
 * @param[in] iter Fixed number of iterations.
 *
 * If the fixed number of iterations is set, the integration algorithm will
 * always performed the given number of iterations, irrespectively of the
 * precision that is reached. This feature is relevant for computing
 * numerical derivates from numerically integrated functions.
 ***************************************************************************/
inline
void GIntegral::fixed_iter(const int& iter)
{
    m_fix_iter = iter;
    return;
}


/***********************************************************************//**
 * @brief Return fixed number of iterations
 *
 * @return Fixed number of iterations.
 ***************************************************************************/
inline
const int& GIntegral::fixed_iter(void) const
{
    return m_fix_iter;
}


/***********************************************************************//**
 * @brief Set relative precision
 *
 * @param[in] eps Relative precision.
 ***************************************************************************/
inline
void GIntegral::eps(const double& eps)
{
    m_eps = eps;
    return;
}


/***********************************************************************//**
 * @brief Get relative precision
 *
 * @return Relative precision.
 ***************************************************************************/
inline
const double& GIntegral::eps(void) const
{
    return m_eps;
}


/***********************************************************************//**
 * @brief Get number of function calls
 *
 * @return Number of function calls.
 ***************************************************************************/
inline
const int& GIntegral::calls(void) const
{
    return m_calls;
}


/***********************************************************************//**
 * @brief Set silence flag
 *
 * @param[in] silent Silence flag.
 ***************************************************************************/
inline
void GIntegral::silent(const bool& silent)
{
    m_silent = silent;
    return;
}


/***********************************************************************//**
 * @brief Get silence flag
 *
 * @return True is class is silent, false otherwise.
 ***************************************************************************/
inline
const bool& GIntegral::silent(void) const
{
    return m_silent;
}


/***********************************************************************//**
 * @brief Set kernel
 *
 * @param[in] kernel Kernel.
 *
 * Sets the kernel for which the integral should be determined.
 ***************************************************************************/
inline
void GIntegral::kernel(GFunction* kernel)
{
    m_kernel = kernel;
    return;
}


/***********************************************************************//**
 * @brief Get kernel
 *
 * @return Kernel.
 ***************************************************************************/
inline
const GFunction* GIntegral::kernel(void) const
{
    return m_kernel;
}


/***********************************************************************//**
 * @brief Signal if integration result is valid
 *
 * @return True is integration result is valid.
 ***************************************************************************/
inline
const bool& GIntegral::is_valid(void) const
{
    return m_isvalid;
}


/***********************************************************************//**
 * @brief Return integration status message
 *
 * @return Integration status message.
 ***************************************************************************/
inline
const std::string& GIntegral::message(void) const
{
    return m_message;
}

#endif /* GINTEGRAL_HPP */
