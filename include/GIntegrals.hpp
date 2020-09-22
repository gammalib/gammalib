/***************************************************************************
 *          GIntegrals.hpp - Integration class for set of functions        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GIntegrals.hpp
 * @brief Integration class for set of functions interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GINTEGRALS_HPP
#define GINTEGRALS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GBase.hpp"
#include "GVector.hpp"

/* __ Forward declarations _______________________________________________ */
class GFunctions;

/***********************************************************************//**
 * @class GIntegrals
 *
 * @brief Integration class for set of functions
 *
 * This class allows to perform integration of a set of functions. The
 * integrand is implemented by a derived class of GFunctions.
 ***************************************************************************/
class GIntegrals : public GBase {

public:

    // Constructors and destructors
    explicit GIntegrals(void);
    explicit GIntegrals(GFunctions* kernels);
    GIntegrals(const GIntegrals& integral);
    virtual ~GIntegrals(void);

    // Operators
    GIntegrals& operator=(const GIntegrals& integral);

    // Methods
    void               clear(void);
    GIntegrals*        clone(void) const;
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
    void               kernels(GFunctions* kernels);
    const GFunctions*  kernels(void) const;
    GVector            romberg(std::vector<double> bounds,
                               const int&          order = 5);
    GVector            romberg(const double& a,
                               const double& b,
                               const int&    order = 5);
    GVector            trapzd(const double& a,
                              const double& b,
                              const int&    n,
                              GVector       result);
    std::string        print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void   init_members(void);
    void   copy_members(const GIntegrals& integral);
    void   free_members(void);
    double polint(const double*  xa,
                  const GVector* ya,
                  const int&     n,
                  const int&     index,
                  const double&  x,
                  double*        dy);

    // Protected data area
    GFunctions* m_kernels;   //!< Pointer to function kernels
    double      m_eps;       //!< Requested relative integration precision
    int         m_max_iter;  //!< Maximum number of iterations
    int         m_fix_iter;  //!< Fixed number of iterations
    bool        m_silent;    //!< Suppress integration warnings in console

    // Integrator results
    mutable int         m_iter;       //!< Number of iterations used
    mutable int         m_calls;      //!< Number of function calls used
    mutable bool        m_isvalid;    //!< Integration result valid (true=yes)
    mutable bool        m_has_abserr; //!< Has absolute integration errors
    mutable bool        m_has_relerr; //!< Has relative integration errors
    mutable GVector     m_abserr;     //!< Absolute integration errors
    mutable GVector     m_relerr;     //!< Absolute integration errors
    mutable std::string m_message;    //!< Status message (if result is invalid)
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GIntegrals").
 ***************************************************************************/
inline
std::string GIntegrals::classname(void) const
{
    return ("GIntegrals");
}


/***********************************************************************//**
 * @brief Return number of iterations
 *
 * @return Number of iterations.
 ***************************************************************************/
inline
const int& GIntegrals::iter(void) const
{
    return m_iter;
}


/***********************************************************************//**
 * @brief Set maximum number of iterations
 *
 * @param[in] iter Maximum number of iterations.
 ***************************************************************************/
inline
void GIntegrals::max_iter(const int& iter)
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
const int& GIntegrals::max_iter(void) const
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
void GIntegrals::fixed_iter(const int& iter)
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
const int& GIntegrals::fixed_iter(void) const
{
    return m_fix_iter;
}


/***********************************************************************//**
 * @brief Set relative precision
 *
 * @param[in] eps Relative precision.
 ***************************************************************************/
inline
void GIntegrals::eps(const double& eps)
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
const double& GIntegrals::eps(void) const
{
    return m_eps;
}


/***********************************************************************//**
 * @brief Get number of function calls
 *
 * @return Number of function calls.
 ***************************************************************************/
inline
const int& GIntegrals::calls(void) const
{
    return m_calls;
}


/***********************************************************************//**
 * @brief Set silence flag
 *
 * @param[in] silent Silence flag.
 ***************************************************************************/
inline
void GIntegrals::silent(const bool& silent)
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
const bool& GIntegrals::silent(void) const
{
    return m_silent;
}


/***********************************************************************//**
 * @brief Set function kernels
 *
 * @param[in] kernels Function kernels.
 *
 * Sets the function kernels for which the integral should be determined.
 ***************************************************************************/
inline
void GIntegrals::kernels(GFunctions* kernels)
{
    m_kernels = kernels;
    return;
}


/***********************************************************************//**
 * @brief Get function kernels
 *
 * @return Function kernels.
 ***************************************************************************/
inline
const GFunctions* GIntegrals::kernels(void) const
{
    return m_kernels;
}


/***********************************************************************//**
 * @brief Signal if integration result is valid
 *
 * @return True is integration result is valid.
 ***************************************************************************/
inline
const bool& GIntegrals::is_valid(void) const
{
    return m_isvalid;
}


/***********************************************************************//**
 * @brief Return integration status message
 *
 * @return Integration status message.
 ***************************************************************************/
inline
const std::string& GIntegrals::message(void) const
{
    return m_message;
}

#endif /* GINTEGRALS_HPP */
