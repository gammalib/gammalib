/***************************************************************************
 *      GOptimizerFunction.i - Optimizer function abstract base class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Juergen Knoedlseder                              *
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
 * @file GOptimizerFunction.i
 * @brief Optimizer function abstract base class
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GOptimizerFunction.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GOptimizerFunction
 *
 * @brief Optimizer function abstract base class
 ***************************************************************************/
class GOptimizerFunction {
public:
    // Constructors and destructors
    GOptimizerFunction(void);
    GOptimizerFunction(const GOptimizerFunction& fct);
    virtual ~GOptimizerFunction(void);

    // Virtual methods
    virtual void           eval(const GOptimizerPars& pars) = 0;
    virtual double         value(void) = 0;
    virtual GVector*       gradient(void) = 0;
    virtual GMatrixSparse* curvature(void) = 0;
};


/***********************************************************************//**
 * @brief GOptimizerFunction class extension
 ***************************************************************************/
%extend GOptimizerFunction {
};
