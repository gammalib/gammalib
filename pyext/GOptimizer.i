/***************************************************************************
 *             GOptimizer.i - Abstract base class for optimizer            *
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
 * @file GOptimizer.i
 * @brief Abstract optimizer abstract base class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GOptimizer.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GOptimizer
 *
 * @brief Abstract optimizer abstract base class
 ***************************************************************************/
class GOptimizer : public GBase {

public:
    // Constructors and destructors
    GOptimizer(void);
    GOptimizer(const GOptimizer& opt);
    virtual ~GOptimizer(void);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GOptimizer* clone(void) const = 0;
    virtual void        optimize(GOptimizerFunction& fct, GOptimizerPars& pars) = 0;
    virtual double      value(void) const = 0;
    virtual int         status(void) const = 0;
    virtual int         iter(void) const = 0;
};


/***********************************************************************//**
 * @brief GOptimizer class extension
 ***************************************************************************/
%extend GOptimizer {
};
