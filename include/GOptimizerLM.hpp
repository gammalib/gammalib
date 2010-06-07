/***************************************************************************
 *            GOptimizerLM.hpp  -  Levenberg Marquardt optimizer           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009 by Jurgen Knodlseder                                *
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
#include "GOptimizer.hpp"
#include "GOptimizerFunction.hpp"
#include "GModels.hpp"


/* __ Includes ___________________________________________________________ */
#define G_LM_CONVERGED            0
#define G_LM_STALLED              1
#define G_LM_SINGULAR             2
#define G_LM_NOT_POSTIVE_DEFINITE 3

/***********************************************************************//**
 * @class GOptimizerLM
 *
 * @brief GOptimizerLM class interface defintion.
 ***************************************************************************/
class GOptimizerLM : public GOptimizer {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GOptimizerLM& opt);

public:

    // Constructors and destructors
    GOptimizerLM();
    GOptimizerLM(const GOptimizerLM& opt);
    virtual ~GOptimizerLM();

    // Operators
    GOptimizerLM&   operator= (const GOptimizerLM& opt);
    GOptimizerPars& operator() (GOptimizerFunction& fct, GOptimizerPars& p);
    GModels&        operator() (GOptimizerFunction& fct, GModels& m);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GOptimizerLM& opt);
    void free_members(void);
    void optimize(GOptimizerFunction* fct, GOptimizerPars* pars);
    void iteration(GOptimizerFunction* fct, GOptimizerPars* pars);

    // Protected data area
    double m_lambda_start;      //!< Initial start value
    double m_lambda_inc;        //!< Lambda increase
    double m_lambda_dec;        //!< Lambda decrease
    double m_eps;               //!< Absolute precision
    int    m_max_iter;          //!< Maximum number of iterations
    int    m_max_stall;         //!< Maximum number of stalls
    double m_lambda;            //!< Actual lambda
    double m_value;             //!< Actual function value
    int    m_status;            //!< Fit status
    int    m_iter;              //!< Iteration

};

#endif /* GOPTIMIZERLM_HPP */
