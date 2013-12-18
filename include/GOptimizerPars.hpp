/***************************************************************************
 *         GOptimizerPars.hpp - Optimizer parameter container class        *
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
 * @brief Optimizer parameters base class definition
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
    explicit GOptimizerPars(const int& number);
    GOptimizerPars(const GOptimizerPars& pars);
    virtual ~GOptimizerPars(void);

    // Operators
    GOptimizerPars&      operator=(const GOptimizerPars& pars);
    GOptimizerPar*       operator[](const int& index);
    const GOptimizerPar* operator[](const int& index) const;
    GOptimizerPar*       operator[](const std::string& name);
    const GOptimizerPar* operator[](const std::string& name) const;

    // Methods
    void            clear(void);
    GOptimizerPars* clone(void) const;
    int             size(void) const;
    bool            isempty(void) const;
    int             nfree(void) const;
    void            attach(GOptimizerPar *par);
    void            remove(const int& index);
    void            reserve(const int& num);
    std::string     print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GOptimizerPars& pars);
    void free_members(void);
    int  get_index(const std::string& name) const;

    // Proteced members
    std::vector<GOptimizerPar*> m_pars;   //!< List of parameters
    std::vector<bool>           m_alloc;  //!< Flags allocation
};


/***********************************************************************//**
 * @brief Return pointer to model
 *
 * @param[in] index Model index [0,...,size()-1].
 *
 * Returns a pointer to the model with the specified @p index.
 ***************************************************************************/
inline
GOptimizerPar* GOptimizerPars::operator[](const int& index)
{
    return (m_pars[index]);
}


/***********************************************************************//**
 * @brief Return pointer to model (const version)
 *
 * @param[in] index Model index [0,...,size()-1].
 *
 * Returns a const pointer to the model with the specified @p index.
 ***************************************************************************/
inline
const GOptimizerPar* GOptimizerPars::operator[](const int& index) const
{
    return (m_pars[index]);
}


/***********************************************************************//**
 * @brief Return number of models in container
 *
 * @return Number of models in container.
 *
 * Returns the number of models in the model container.
 ***************************************************************************/
inline
int GOptimizerPars::size(void) const
{
    return (m_pars.size());
}


/***********************************************************************//**
 * @brief Signals if there are no models in container
 *
 * @return True if container is empty, false otherwise.
 *
 * Signals if the model container does not contain any model.
 ***************************************************************************/
inline
bool GOptimizerPars::isempty(void) const
{
    return (m_pars.empty());
}


/***********************************************************************//**
 * @brief Reserves space for models in container
 *
 * @param[in] num Number of models
 *
 * Reserves space for @p num models in the container.
 ***************************************************************************/
inline
void GOptimizerPars::reserve(const int& num)
{
    m_pars.reserve(num);
    return;
}

#endif /* GOPTIMIZERPARS_HPP */
