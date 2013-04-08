/***************************************************************************
 *             GOptimizerLM.hpp - Levenberg Marquardt optimizer            *
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
    explicit GOptimizerLM(GLog& log);
    GOptimizerLM(const GOptimizerLM& opt);
    virtual ~GOptimizerLM(void);

    // Operators
    GOptimizerLM& operator=(const GOptimizerLM& opt);

    // Implemented pure virtual base class methods
    virtual void          clear(void);
    virtual GOptimizerLM* clone(void) const;
    virtual void          optimize(GOptimizerFunction& fct, GOptimizerPars& pars);
    virtual double        value(void) const { return m_value; }   //!< @brief Return function value
    virtual int           status(void) const { return m_status; } //!< @brief Return optimization status
    virtual int           iter(void) const { return m_iter; }     //!< @brief Return number of iterations
    virtual std::string   print(const GChatter& chatter = NORMAL) const;
    
    // Methods
    void          max_iter(const int& n) { m_max_iter=n; }                //!< @brief Set maximum number of iterations
    void          max_stalls(const int& n) { m_max_stall=n; }             //!< @brief Set maximum number of stalls
    void          max_boundary_hits(const int& n) { m_max_stall=n; }      //!< @brief Set maximum number of boundary hits
    void          lambda_start(const double& val) { m_lambda_start=val; } //!< @brief Set lambda start value
    void          lambda_inc(const double& val) { m_lambda_inc=val; }     //!< @brief Set lambda increment
    void          lambda_dec(const double& val) { m_lambda_dec=val; }     //!< @brief Set lambda decrement
    void          eps(const double& eps) { m_eps=eps; }                   //!< @brief Set convergence precisions
    int           max_iter(void) const { return m_max_iter; }             //!< @brief Return maximum number of iterations
    int           max_stalls(void) const { return m_max_stall; }          //!< @brief Return maximum number of stalls 
    int           max_boundary_hits(void) const { return m_max_hit; }     //!< @brief Return maximum number of boundary hits
    const double& lambda_start(void) const { return m_lambda_start; }     //!< @brief Return lambda start value
    const double& lambda_inc(void) const { return m_lambda_inc; }         //!< @brief Return lambda increment
    const double& lambda_dec(void) const { return m_lambda_dec; }         //!< @brief Return lambda derement
    const double& lambda(void) const { return m_lambda; }                 //!< @brief Return lambda value
    const double& eps(void) const { return m_eps; }                       //!< @brief Return convergence precision

protected:
    // Protected methods
    void   init_members(void);
    void   copy_members(const GOptimizerLM& opt);
    void   free_members(void);
    void   iteration(GOptimizerFunction& fct, GOptimizerPars& pars);
    void   errors(GOptimizerFunction& fct, GOptimizerPars& pars);
    double step_size(const GVector& grad, const GOptimizerPars& pars);

    // Protected members
    int               m_npars;           //!< Number of parameters
    int               m_nfree;           //!< Number of free parameters
    double            m_lambda_start;    //!< Initial start value
    double            m_lambda_inc;      //!< Lambda increase
    double            m_lambda_dec;      //!< Lambda decrease
    double            m_eps;             //!< Absolute precision
    int               m_max_iter;        //!< Maximum number of iterations
    int               m_max_stall;       //!< Maximum number of stalls
    int               m_max_hit;         //!< Maximum number of successive hits
    bool              m_step_adjust;     //!< Adjust step size to boundaries
    std::vector<bool> m_hit_boundary;    //!< Bookkeeping array for boundary hits
    std::vector<int>  m_hit_minimum;     //!< Bookkeeping of successive minimum hits
    std::vector<int>  m_hit_maximum;     //!< Bookkeeping of successive maximum hits
    std::vector<bool> m_par_freeze;      //!< Bookkeeping of parameter freeze
    std::vector<bool> m_par_remove;      //!< Bookkeeping of parameter removal
    double            m_lambda;          //!< Actual lambda
    double            m_value;           //!< Actual function value
    int               m_status;          //!< Fit status
    int               m_iter;            //!< Iteration
    GLog*             m_logger;          //!< Pointer to optional logger

};

#endif /* GOPTIMIZERLM_HPP */
