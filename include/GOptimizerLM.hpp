/***************************************************************************
 *             GOptimizerLM.hpp - Levenberg Marquardt optimizer            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2017 by Juergen Knoedlseder                         *
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
 * @file GOptimizerLM.hpp
 * @brief Levenberg Marquardt optimizer class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GOPTIMIZERLM_HPP
#define GOPTIMIZERLM_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include "GOptimizer.hpp"
#include "GOptimizerFunction.hpp"
#include "GLog.hpp"

/* __ Definitions ________________________________________________________ */
#define G_LM_CONVERGED            0
#define G_LM_STALLED              1
#define G_LM_SINGULAR             2
#define G_LM_NOT_POSTIVE_DEFINITE 3
#define G_LM_BAD_ERRORS           4


/***********************************************************************//**
 * @class GOptimizerLM
 *
 * @brief Levenberg Marquardt optimizer class
 *
 * This method implements an Levenberg Marquardt optimizer.
 ***************************************************************************/
class GOptimizerLM : public GOptimizer {

public:

    // Constructors and destructors
    GOptimizerLM(void);
    explicit GOptimizerLM(GLog* log);
    GOptimizerLM(const GOptimizerLM& opt);
    virtual ~GOptimizerLM(void);

    // Operators
    GOptimizerLM& operator=(const GOptimizerLM& opt);

    // Implemented pure virtual base class methods
    virtual void          clear(void);
    virtual GOptimizerLM* clone(void) const;
    virtual std::string   classname(void) const;
    virtual void          optimize(GOptimizerFunction& fct, GOptimizerPars& pars);
    virtual void          errors(GOptimizerFunction& fct, GOptimizerPars& pars);
    virtual double        value(void) const;
    virtual int           status(void) const;
    virtual int           iter(void) const;
    virtual std::string   print(const GChatter& chatter = NORMAL) const;
    
    // Methods
    void          logger(GLog* log);
    void          max_iter(const int& max_iter);
    void          max_stalls(const int& max_stalls);
    void          max_boundary_hits(const int& max_hit);
    void          lambda_start(const double& value);
    void          lambda_inc(const double& value);
    void          lambda_dec(const double& value);
    void          eps(const double& eps);
    void          accept_dec(const double& value);
    int           max_iter(void) const;
    int           max_stalls(void) const;
    int           max_boundary_hits(void) const;
    const double& lambda_start(void) const;
    const double& lambda_inc(void) const;
    const double& lambda_dec(void) const;
    const double& lambda(void) const;
    const double& eps(void) const;
    const double& accept_dec(void) const;

protected:
    // Protected methods
    void   init_members(void);
    void   copy_members(const GOptimizerLM& opt);
    void   free_members(void);
    double iteration(GOptimizerFunction& fct, GOptimizerPars& pars);
    double step_size(const GVector& grad, const GOptimizerPars& pars);

    // Protected members
    int               m_npars;        //!< Number of parameters
    int               m_nfree;        //!< Number of free parameters
    double            m_lambda_start; //!< Initial start value
    double            m_lambda_inc;   //!< Lambda increase
    double            m_lambda_dec;   //!< Lambda decrease
    double            m_eps;          //!< Absolute precision
    double            m_accept_dec;   //!< Acceptable function decrease
    int               m_max_iter;     //!< Maximum number of iterations
    int               m_max_stall;    //!< Maximum number of stalls
    int               m_max_hit;      //!< Maximum number of successive hits
    int               m_max_dec;      //!< Maximum number of function decrease
    bool              m_step_adjust;  //!< Adjust step size to boundaries
    std::vector<bool> m_hit_boundary; //!< Bookkeeping array for boundary hits
    std::vector<int>  m_hit_minimum;  //!< Bookkeeping of successive minimum hits
    std::vector<int>  m_hit_maximum;  //!< Bookkeeping of successive maximum hits
    std::vector<bool> m_par_freeze;   //!< Bookkeeping of parameter freeze
    std::vector<bool> m_par_remove;   //!< Bookkeeping of parameter removal
    double            m_lambda;       //!< Actual lambda
    double            m_value;        //!< Actual function value
    double            m_delta;        //!< Function improvement
    int               m_status;       //!< Fit status
    int               m_iter;         //!< Iteration
    int               m_num_dec;      //!< Number of function decreases
    GLog*             m_logger;       //!< Pointer to optional logger

};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GOptimizerLM").
 ***************************************************************************/
inline
std::string GOptimizerLM::classname(void) const
{
    return ("GOptimizerLM");
}


/***********************************************************************//**
 * @brief Return function value
 *
 * @return Function value.
 ***************************************************************************/
inline
double GOptimizerLM::value(void) const
{
    return (m_value);
}


/***********************************************************************//**
 * @brief Return optimizer status
 *
 * @return Optimizer status.
 ***************************************************************************/
inline
int GOptimizerLM::status(void) const
{
    return (m_status);
}


/***********************************************************************//**
 * @brief Return number of iterations
 *
 * @return Number of iterations.
 ***************************************************************************/
inline
int GOptimizerLM::iter(void) const
{
    return (m_iter);
}


/***********************************************************************//**
 * @brief Set logger
 *
 * @param[in] log Logger to use in optimizer.
 *
 * Set the logger into which the optimizer will dump any output.
 ***************************************************************************/
inline
void GOptimizerLM::logger(GLog* log)
{
    m_logger = log;
    return;
}


/***********************************************************************//**
 * @brief Set maximum number of iterations
 *
 * @param[in] max_iter Maximum number of iterations.
 ***************************************************************************/
inline
void GOptimizerLM::max_iter(const int& max_iter)
{
    m_max_iter = max_iter;
    return;
}


/***********************************************************************//**
 * @brief Set maximum number of allowed subsequent stalls
 *
 * @param[in] max_stalls Maximum number of allowed subsequent stalls.
 ***************************************************************************/
inline
void GOptimizerLM::max_stalls(const int& max_stalls)
{
    m_max_stall = max_stalls;
    return;
}


/***********************************************************************//**
 * @brief Set maximum number of parameter boundary hits
 *
 * @param[in] max_hit Maximum number of parameter boundary hits.
 ***************************************************************************/
inline
void GOptimizerLM::max_boundary_hits(const int& max_hit)
{
    m_max_hit = max_hit;
    return;
}


/***********************************************************************//**
 * @brief Set lambda starting value
 *
 * @param[in] value Lambda starting value.
 ***************************************************************************/
inline
void GOptimizerLM::lambda_start(const double& value)
{
    m_lambda_start = value;
    return;
}


/***********************************************************************//**
 * @brief Set lambda increment value
 *
 * @param[in] value Lambda increment value.
 ***************************************************************************/
inline
void GOptimizerLM::lambda_inc(const double& value)
{
    m_lambda_inc = value;
    return;
}


/***********************************************************************//**
 * @brief Set lambda decrement value
 *
 * @param[in] value Lambda decrement value.
 ***************************************************************************/
inline
void GOptimizerLM::lambda_dec(const double& value)
{
    m_lambda_dec = value;
    return;
}


/***********************************************************************//**
 * @brief Set requested absolute convergence precision
 *
 * @param[in] eps Requested absolute convergence precision.
 ***************************************************************************/
inline
void GOptimizerLM::eps(const double& eps)
{
    m_eps = eps;
    return;
}


/***********************************************************************//**
 * @brief Set acceptable function decrease
 *
 * @param[in] value Acceptable function decrease.
 *
 * Sets the acceptable function decrease value for which the new solution
 * will be kept and the iterations continue. This strategy provides better
 * convergence in case that a function decrease is encountered.
 ***************************************************************************/
inline
void GOptimizerLM::accept_dec(const double& value)
{
    m_accept_dec = value;
    return;
}




/***********************************************************************//**
 * @brief Return maximum number of iterations
 *
 * @return Maximum number of iterations.
 ***************************************************************************/
inline
int GOptimizerLM::max_iter(void) const
{
    return (m_max_iter);
}


/***********************************************************************//**
 * @brief Return maximum number of allowed subsequent stalls
 *
 * @return Maximum number of allowed subsequent stalls.
 ***************************************************************************/
inline
int GOptimizerLM::max_stalls(void) const
{
    return (m_max_stall);
}


/***********************************************************************//**
 * @brief Return maximum number of parameter boundary hits
 *
 * @return Maximum number of parameter boundary hits.
 ***************************************************************************/
inline
int GOptimizerLM::max_boundary_hits(void) const
{
    return (m_max_hit);
}


/***********************************************************************//**
 * @brief Return lambda starting value
 *
 * @return Lambda starting value.
 ***************************************************************************/
inline
const double& GOptimizerLM::lambda_start(void) const
{
    return (m_lambda_start);
}


/***********************************************************************//**
 * @brief Return lambda increment value
 *
 * @return Lambda increment value.
 ***************************************************************************/
inline
const double& GOptimizerLM::lambda_inc(void) const
{
    return (m_lambda_inc);
}


/***********************************************************************//**
 * @brief Return lambda decrement value
 *
 * @return Lambda decrement value.
 ***************************************************************************/
inline
const double& GOptimizerLM::lambda_dec(void) const
{
    return (m_lambda_dec);
}


/***********************************************************************//**
 * @brief Return lambda value
 *
 * @return Lambda value.
 ***************************************************************************/
inline
const double& GOptimizerLM::lambda(void) const
{
    return (m_lambda);
}


/***********************************************************************//**
 * @brief Return requested absolute convergence precision
 *
 * @return Requested absolute convergence precision.
 ***************************************************************************/
inline
const double& GOptimizerLM::eps(void) const
{
    return (m_eps);
}


/***********************************************************************//**
 * @brief Return acceptable function decrease
 *
 * @return Acceptable function decrease.
 ***************************************************************************/
inline
const double& GOptimizerLM::accept_dec(void) const
{
    return (m_accept_dec);
}

#endif /* GOPTIMIZERLM_HPP */
