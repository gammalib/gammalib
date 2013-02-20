/***************************************************************************
 *    GOptimizerPars.hpp - Abstract optimizer parameters container class   *
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
 * @file GOptimizerPars.hpp
 * @brief Abstract optimizer parameters base class definition
 * @author Juergen Knoedlseder
 */

#ifndef GOPTIMIZERPARS_HPP
#define GOPTIMIZERPARS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GContainer.hpp"
#include "GModelPar.hpp"


/***********************************************************************//**
 * @class GOptimizerPars
 *
 * @brief Optimizer parameter container class
 *
 * This class is a container class for parameters of a function that are to
 * be optimized. The optimizer function is defined by the abstract
 * GOptimizerFunction class.
 *
 * The class holds a flat array of pointers to models parameters. It neither
 * deals with allocation and deallocation, nor with cloning of function
 * parameters, this will be done by the classes that actually implement the
 * model parameters.
 *
 * @todo This container class has no operator[] method as GModels is a
 *       derived class of this container class so that GModels can be
 *       passed to the optimizer. GModels has an operator[] to access
 *       models, so we cannot implement this here again. This is not
 *       very clean. We should think about how we can make this a clean
 *       thing ...
 ***************************************************************************/
class GOptimizerPars : public GContainer {

public:
    // Constructors and destructors
    GOptimizerPars(void);
    GOptimizerPars(const GOptimizerPars& pars);
    virtual ~GOptimizerPars(void);

    // Operators
    virtual GOptimizerPars& operator=(const GOptimizerPars& pars);

    // Pure virtual base class methods
    virtual void             clear(void) = 0;
    virtual GOptimizerPars*  clone(void) const = 0;
    virtual int              size(void) const = 0;
    virtual bool             isempty(void) const = 0;
    virtual void             pop(const int& index) = 0;
    virtual void             reserve(const int& num) = 0;
    virtual std::string      print(void) const = 0;

    // Other methods
    virtual int              npars(void) const { return m_pars.size(); } //! @brief Return number of parameters
    virtual int              nfree(void) const;
    virtual GModelPar&       par(const int& index);
    virtual const GModelPar& par(const int& index) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GOptimizerPars& pars);
    void free_members(void);

    // Proteced members
    std::vector<GModelPar*> m_pars;   //!< Pointers to model parameters
};

#endif /* GOPTIMIZERPARS_HPP */
