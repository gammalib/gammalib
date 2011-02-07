/***************************************************************************
 *            GOptimizerLM.hpp  -  Levenberg Marquardt optimizer           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GOptimizerLM.hpp
 * @brief GOptimizerLM base class interface definition.
 * @author J. Knodlseder
 */

#ifndef GOPTIMIZERLM_HPP
#define GOPTIMIZERLM_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include "GOptimizer.hpp"
#include "GOptimizerFunction.hpp"
#include "GModels.hpp"
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
 * @brief GOptimizerLM class interface defintion.
 ***************************************************************************/
class GOptimizerLM : public GOptimizer {

public:

    // Constructors and destructors
    GOptimizerLM(void);
    GOptimizerLM(GLog& log);
    GOptimizerLM(const GOptimizerLM& opt);
    virtual ~GOptimizerLM(void);

    // Operators
    GOptimizerLM&   operator= (const GOptimizerLM& opt);
    GOptimizerPars& operator() (GOptimizerFunction& fct, GOptimizerPars& p);
    GModels&        operator() (GOptimizerFunction& fct, GModels& m);

    // Implemented pure virtual methods
    void          clear(void);
    GOptimizerLM* clone(void) const;
    std::string   print(void) const;
    
    // Methods
    void        max_iter(const int& n) { m_max_iter=n; }
    void        max_stalls(const int& n) { m_max_stall=n; }
    void        max_boundary_hits(const int& n) { m_max_stall=n; }
    void        lambda_start(const double& val) { m_lambda_start=val; }
    void        lambda_inc(const double& val) { m_lambda_inc=val; }
    void        lambda_dec(const double& val) { m_lambda_dec=val; }
    void        eps(const double& eps) { m_eps=eps; }
    int         max_iter(void) const { return m_max_iter; }
    int         max_stalls(void) const { return m_max_stall; }
    int         max_boundary_hits(void) const { return m_max_hit; }
    int         status(void) const { return m_status; }
    int         iter(void) const { return m_iter; }
    double      lambda_start(void) const { return m_lambda_start; }
    double      lambda_inc(void) const { return m_lambda_inc; }
    double      lambda_dec(void) const { return m_lambda_dec; }
    double      lambda(void) const { return m_lambda; }
    double      eps(void) const { return m_eps; }
    double      value(void) const { return m_value; }

protected:
    // Protected methods
    void   init_members(void);
    void   copy_members(const GOptimizerLM& opt);
    void   free_members(void);
    void   optimize(GOptimizerFunction* fct, GOptimizerPars* pars);
    void   iteration(GOptimizerFunction* fct, GOptimizerPars* pars);
    void   errors(GOptimizerFunction* fct, GOptimizerPars* pars);
    double step_size(GVector* grad, GOptimizerPars* pars);

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
    double            m_lambda;          //!< Actual lambda
    double            m_value;           //!< Actual function value
    int               m_status;          //!< Fit status
    int               m_iter;            //!< Iteration
    GLog*             m_logger;          //!< Pointer to optional logger

};

#endif /* GOPTIMIZERLM_HPP */
