/***************************************************************************
 *                   GDerivative.hpp - Derivative class                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2014 by Juergen Knoedlseder                         *
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
 * @file GDerivative.hpp
 * @brief GDerivative class interface definition.
 * @author Juergen Knoedlseder
 */

#ifndef GDERIVATIVE_HPP
#define GDERIVATIVE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GFunction.hpp"


/***********************************************************************//**
 * @class GDerivative
 *
 * @brief Numerical derivatives class
 *
 * This class allows to compute numerical derivatives using various methods.
 * The function to be derived is implemented by the abstract GFunction
 * class.
 ***************************************************************************/
class GDerivative : public GBase {

public:

    // Constructors and destructors
    GDerivative(void);
    explicit GDerivative(GFunction* func);
    GDerivative(const GDerivative& dx);
    virtual ~GDerivative(void);

    // Operators
    GDerivative& operator=(const GDerivative& dx);

    // Methods
    void             clear(void);
    GDerivative*     clone(void) const;
    void             max_iter(const int& max_iter);
    void             eps(const double& eps);
    void             step_frac(const double& fraction);
    void             silent(const bool& silent);
    const int&       iter(void) const;
    const int&       max_iter(void) const;
    const double&    eps(void) const;
    const double&    step_frac(void) const;
    const bool&      silent(void) const;
    void             function(GFunction* func);
    const GFunction* function(void) const;
    double           value(const double& x, const double& step = 0.0);
    double           ridder(const double& x, const double& h, double* err);
    double           minuit2(const double& x, double* err);
    double           difference(const double& x, const double& h);
    double           left_difference(const double& x, const double& h);
    double           right_difference(const double& x, const double& h);
    double           smooth_robust(const double& x, const double& h, 
                                   const int& degree = 2,
                                   const int& length = 5);
    std::string      print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GDerivative& dx);
    void free_members(void);
    void set_tiny(void);

    // Tiny number computation
    class tiny {
    public:
        tiny(void) : m_one(1.0) {}
        ~tiny(void) {}
        double one(void) const;
        double operator()(double eps) const;
    private:
        double m_one;
    };

    // Protected members
    GFunction* m_func;         //!< Pointer to function
    double     m_eps;          //!< Derivative precision
    double     m_step_frac;    //!< Value fraction to use for initial step
    double     m_tiny;         //!< Tiny number for minuit2
    int        m_max_iter;     //!< Maximum number of iterations
    int        m_iter;         //!< Number of iterations used
    bool       m_silent;       //!< Suppress warnings
};


/***********************************************************************//**
 * @brief Return number of iterations
 *
 * @return Number of iterations.
 ***************************************************************************/
inline
const int& GDerivative::iter(void) const
{
    return m_iter;
}


/***********************************************************************//**
 * @brief Set maximum number of iterations
 *
 * @param[in] max_iter Maximum number of iterations.
 ***************************************************************************/
inline
void GDerivative::max_iter(const int& max_iter)
{
    m_max_iter = max_iter;
    return;
}


/***********************************************************************//**
 * @brief Return maximum number of iterations
 *
 * @return Maximum number of iterations.
 ***************************************************************************/
inline
const int& GDerivative::max_iter(void) const
{
    return m_max_iter;
}


/***********************************************************************//**
 * @brief Set precision
 *
 * @param[in] eps Precision.
 ***************************************************************************/
inline
void GDerivative::eps(const double& eps)
{
    m_eps = eps;
    return;
}


/***********************************************************************//**
 * @brief Get precision
 *
 * @return Precision.
 ***************************************************************************/
inline
const double& GDerivative::eps(void) const
{
    return m_eps;
}


/***********************************************************************//**
 * @brief Set step fraction
 *
 * @param[in] fraction Step fraction.
 ***************************************************************************/
inline
void GDerivative::step_frac(const double& fraction)
{
    m_step_frac = fraction;
    return;
}


/***********************************************************************//**
 * @brief Get step fraction
 *
 * @return Step fraction.
 ***************************************************************************/
inline
const double& GDerivative::step_frac(void) const
{
    return m_step_frac;
}


/***********************************************************************//**
 * @brief Set silence flag
 *
 * @param[in] silent Silence flag.
 ***************************************************************************/
inline
void GDerivative::silent(const bool& silent)
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
const bool& GDerivative::silent(void) const
{
    return m_silent;
}


/***********************************************************************//**
 * @brief Set function
 *
 * @param[in] function Function.
 *
 * Sets the function for which the derivative should be determined.
 ***************************************************************************/
inline
void GDerivative::function(GFunction* function)
{
    m_func = function;
    return;
}


/***********************************************************************//**
 * @brief Get function
 *
 * @return Function.
 ***************************************************************************/
inline
const GFunction* GDerivative::function(void) const
{
    return m_func;
}

#endif /* GDERIVATIVE_HPP */
