/***************************************************************************
 * GOptimizerLM.i  -  Levenberg Marquardt optimizer class Python interface *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2012 by Juergen Knoedlseder                         *
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
 * @file GOptimizerLM.i
 * @brief Levenberg Marquardt optimizer class Pyhton interface definition
 * @author Juergen Knoedlseder
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

    // Implemented pure virtual methods
    virtual void          clear(void);
    virtual GOptimizerLM* clone(void) const;
    virtual void          optimize(GOptimizerFunction& fct, GOptimizerPars& pars);
    virtual double        value(void) const;
    virtual int           status(void) const;
    virtual int           iter(void) const;
    
    // Methods
    void          max_iter(const int& n);
    void          max_stalls(const int& n);
    void          max_boundary_hits(const int& n);
    void          lambda_start(const double& val);
    void          lambda_inc(const double& val);
    void          lambda_dec(const double& val);
    void          eps(const double& eps);
    int           max_iter(void) const;
    int           max_stalls(void) const;
    int           max_boundary_hits(void) const;
    const double& lambda_start(void) const;
    const double& lambda_inc(void) const;
    const double& lambda_dec(void) const;
    const double& eps(void) const;
};


/***********************************************************************//**
 * @brief GOptimizerLM class extension
 ***************************************************************************/
%extend GOptimizerLM {
    GOptimizerLM copy() {
        return (*self);
    }
    const double& lambda_value(void) {
        return (self->lambda());
    }
};
