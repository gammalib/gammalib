/***************************************************************************
 *           GOptimizer.hpp  -  Abstract base class for optimizer          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2011 by Jurgen Knodlseder                           *
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
 * @brief GOptimizer abstract base class interface definition.
 * @author J. Knodlseder
 */

#ifndef GOPTIMIZER_HPP
#define GOPTIMIZER_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GLog.hpp"
#include "GOptimizerPars.hpp"
#include "GOptimizerFunction.hpp"
#include "GModels.hpp"


/***********************************************************************//**
 * @class GOptimizer
 *
 * @brief GOptimizer abstract base class interface defintion.
 ***************************************************************************/
class GOptimizer {

    // I/O friends
    friend std::ostream& operator<<(std::ostream& os, const GOptimizer& opt);
    friend GLog&         operator<<(GLog& log,        const GOptimizer& opt);

public:
    // Constructors and destructors
    GOptimizer(void);
    GOptimizer(const GOptimizer& opt);
    virtual ~GOptimizer();

    // Operators
    virtual GOptimizer&     operator= (const GOptimizer& opt);
    virtual GOptimizerPars& operator() (GOptimizerFunction& fct, GOptimizerPars& p) = 0;
    virtual GModels&        operator() (GOptimizerFunction& fct, GModels& m) = 0;

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GOptimizer* clone(void) const = 0;
    virtual std::string print(void) const = 0;
 
    // Implemented methods

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GOptimizer& opt);
    void free_members(void);
};

#endif /* GOPTIMIZER_HPP */
