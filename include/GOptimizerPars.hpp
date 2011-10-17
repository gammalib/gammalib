/***************************************************************************
 *        GOptimizerPars.hpp  -  Optimizer parameter container class       *
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
 * @file GOptimizerPars.hpp
 * @brief Optimizer parameter container class interface definition
 * @author J. Knodlseder
 */

#ifndef GOPTIMIZERPARS_HPP
#define GOPTIMIZERPARS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GLog.hpp"
#include "GModelPar.hpp"


/***********************************************************************//**
 * @class GOptimizerPars
 *
 * @brief Optimizer parameter container class
 *
 * This class is a container class for model parameters.
 *
 * As the class holds simpliy a collection of model parameters, it should
 * neither deal with allocation and deallocation, nor with cloning of
 * model parameters. This will be done by the classes that actually
 * implement the model parameters.
 *
 * @todo This container class has no operator[] method as GModels is a
 *       derived class of this container class so that GModels can be
 *       passed to the optimizer. GModels has an operator[] to access
 *       models, so we cannot implement this here again. This is not
 *       very clean. We should think about how we can make this a clean
 *       thing ...
 ***************************************************************************/
class GOptimizerPars {

    // I/O friends
    friend std::ostream& operator<<(std::ostream& os, const GOptimizerPars& pars);
    friend GLog&         operator<<(GLog& log,        const GOptimizerPars& pars);

public:
    // Constructors and destructors
    GOptimizerPars(void);
    GOptimizerPars(const GOptimizerPars& pars);
    virtual ~GOptimizerPars(void);

    // Operators
    virtual GOptimizerPars& operator=(const GOptimizerPars& pars);

    // Methods
    int              npars(void) const { return m_pars.size(); }
    int              nfree(void) const;
    GModelPar&       par(const int& index);
    const GModelPar& par(const int& index) const;
    std::string      print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GOptimizerPars& pars);
    void free_members(void);

    // Proteced members
    std::vector<GModelPar*> m_pars;   //!< Pointers to model parameters
};

#endif /* GOPTIMIZERPARS_HPP */
