/***************************************************************************
 *           GOptimizerLM.i - Levenberg Marquardt optimizer class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2020 by Juergen Knoedlseder                         *
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
%}


/***********************************************************************//**
 * @class GOptimizerLM
 *
 * @brief GOptimizerLM class SWIG interface definition.
 ***************************************************************************/
class GOptimizerLM : public GOptimizer {
public:

    // Constructors and destructors
    GOptimizerLM(void);
    explicit GOptimizerLM(GLog* log);
    GOptimizerLM(const GOptimizerLM& opt);
    virtual ~GOptimizerLM(void);

    // Implemented pure virtual methods
    virtual void          clear(void);
    virtual GOptimizerLM* clone(void) const;
    virtual std::string   classname(void) const;
    virtual void          optimize(GOptimizerFunction& fct, GOptimizerPars& pars);
    virtual void          errors(GOptimizerFunction& fct, GOptimizerPars& pars);
    virtual double        value(void) const;
    virtual int           status(void) const;
    virtual int           iter(void) const;
    
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
    const int&    npars(void) const;
    const int&    nfree(void) const;
    const int&    max_iter(void) const;
    const int&    max_stalls(void) const;
    const int&    max_boundary_hits(void) const;
    const double& lambda_start(void) const;
    const double& lambda_inc(void) const;
    const double& lambda_dec(void) const;
    const double& eps(void) const;
    const double& accept_dec(void);
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
%pythoncode {
    def __getstate__(self):
        state = (self.max_iter(), self.max_stalls(), self.max_boundary_hits(),
                 self.lambda_start(), self.lambda_inc(), self.lambda_dec(),
                 self.eps(), self.accept_dec())
        return state
    def __setstate__(self, state):
        self.__init__()
        self.max_iter(state[0])
        self.max_stalls(state[1])
        self.max_boundary_hits(state[2])
        self.lambda_start(state[3])
        self.lambda_inc(state[4])
        self.lambda_dec(state[5])
        self.eps(state[6])
        self.accept_dec(state[7])
}
};
