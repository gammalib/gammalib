/***************************************************************************
 *                      GDerivative.i - Derivative class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
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
 * @file GDerivative.i
 * @brief Derivative class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GDerivative.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GDerivative
 *
 * @brief Numerical derivatives class
 ***************************************************************************/
class GDerivative : public GBase {
public:

    GDerivative(void);
    explicit GDerivative(GFunction* function);
    GDerivative(const GDerivative& dx);
    virtual ~GDerivative(void);

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
    void             function(GFunction* function);
    const GFunction* function(void) const;
    double           value(const double& x, const double& step = 0.0);
    double           ridder(const double& x, const double& h, double* err);
    double           minuit2(const double& x, double* err);
    double           difference(const double& x, const double& h);
};


/***********************************************************************//**
 * @brief GDerivative class extension
 ***************************************************************************/
%extend GDerivative {
    GDerivative copy() {
        return (*self);
    }
};
