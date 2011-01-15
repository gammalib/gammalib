/***************************************************************************
 *  GOptimizerLM.i  -  Levenberg Marquardt optimizer class SWIG interface  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GOptimizerLM.i
 * @brief GOptimizerLM class SWIG interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GOptimizerLM.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GOptimizerLM
 *
 * @brief GOptimizerLM class SWIG interface defintion.
 ***************************************************************************/
class GOptimizerLM : public GOptimizer {
public:

    // Constructors and destructors
    GOptimizerLM(void);
    GOptimizerLM(GLog& log);
    GOptimizerLM(const GOptimizerLM& opt);
    virtual ~GOptimizerLM(void);

    // Operators (we need those to make the class not abstract!!!)
    GOptimizerPars& operator() (GOptimizerFunction& fct, GOptimizerPars& p);
    GModels&        operator() (GOptimizerFunction& fct, GModels& m);

    // Methods
    void   max_iter(const int& n) { m_max_iter=n; return; }
    void   max_stalls(const int& n) { m_max_stall=n; return; }
    void   max_boundary_hits(const int& n) { m_max_stall=n; }
    void   lambda_start(const double& val) { m_lambda_start=val; return; }
    void   lambda_inc(const double& val) { m_lambda_inc=val; return; }
    void   lambda_dec(const double& val) { m_lambda_dec=val; return; }
    void   eps(const double& eps) { m_eps=eps; return; }
    int    max_iter(void) const { return m_max_iter; }
    int    max_stalls(void) const { return m_max_stall; }
    int    max_boundary_hits(void) const { return m_max_hit; }
    int    status(void) const { return m_status; }
    int    iter(void) const { return m_iter; }
    double lambda_start(void) const { return m_lambda_start; }
    double lambda_inc(void) const { return m_lambda_inc; }
    double lambda_dec(void) const { return m_lambda_dec; }
    //double lambda(void) const { return m_lambda; }
    double eps(void) const { return m_eps; }
    double value(void) const { return m_value; }
};


/***********************************************************************//**
 * @brief GOptimizerLM class extension
 ***************************************************************************/
%extend GOptimizerLM {
    char *__str__() {
        return tochar(self->print());
    }
    GOptimizerLM copy() {
        return (*self);
    }
};
