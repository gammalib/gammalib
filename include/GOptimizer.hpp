/***************************************************************************
 *            GOptimizer.hpp - Abstract base class for optimizer           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2014 by Juergen Knoedlseder                         *
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
 * @file GOptimizer.hpp
 * @brief Abstract optimizer abstract base class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GOPTIMIZER_HPP
#define GOPTIMIZER_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GOptimizerPars.hpp"
#include "GOptimizerFunction.hpp"


/***********************************************************************//**
 * @class GOptimizer
 *
 * @brief Abstract optimizer abstract base class
 *
 * This class defines the abstract interface for the optimizer class. The
 * optimizer class is used to optimize the parameters of a function. The
 * function is implemented using the abstract GOptimizerFunction class
 * while the parameters are implemented using the abstract GOptimizerPars
 * class.
 *
 * The main driver of this class is the optimize() method which optimizes
 * the parameters using a function. The function value can be accessed
 * using the value() method, the status() method provides an integer with
 * status information, the iter() method gives the number of iterations for
 * iterative optimization algorithms.
 ***************************************************************************/
class GOptimizer : public GBase {

public:
    // Constructors and destructors
    GOptimizer(void);
    GOptimizer(const GOptimizer& opt);
    virtual ~GOptimizer(void);

    // Operators
    virtual GOptimizer& operator=(const GOptimizer& opt);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GOptimizer* clone(void) const = 0;
    virtual std::string classname(void) const = 0;
    virtual void        optimize(GOptimizerFunction& fct, GOptimizerPars& pars) = 0;
    virtual void        errors(GOptimizerFunction& fct, GOptimizerPars& pars) = 0;
    virtual double      value(void) const = 0;
    virtual int         status(void) const = 0;
    virtual int         iter(void) const = 0;
    virtual std::string print(const GChatter& chatter = NORMAL) const = 0;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GOptimizer& opt);
    void free_members(void);
};

#endif /* GOPTIMIZER_HPP */
