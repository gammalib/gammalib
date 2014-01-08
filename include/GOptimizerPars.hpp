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
#include <vector>
#include "GContainer.hpp"
#include "GOptimizerPar.hpp"


/***********************************************************************//**
 * @class GOptimizerPars
 *
 * @brief Optimizer parameter container class
 *
 * This class is a container class for parameters of a function that are to
 * be optimized. The optimizer function is defined by the abstract
 * GOptimizerFunction class.
 *
 * The class holds a flat array of pointers to optimizer parameters. The
 * class can deal with parameters that are allocated by the class itself,
 * or with pointers that are passed by a client and allocated/deallocated
 * by the client.
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
    void                 clear(void);
    GOptimizerPars*      clone(void) const;
    GOptimizerPar*       at(const int& index);
    const GOptimizerPar* at(const int& index) const;
    int                  size(void) const;
    bool                 is_empty(void) const;
    int                  nfree(void) const;
    GOptimizerPar*       set(const int& index, const GOptimizerPar& par);
    GOptimizerPar*       set(const std::string& name, const GOptimizerPar& par);
    void                 attach(GOptimizerPar* par);
    void                 attach(const int& index, GOptimizerPar* par);
    void                 attach(const std::string& name, GOptimizerPar* par);
    GOptimizerPar*       insert(const int& index, const GOptimizerPar& par);
    GOptimizerPar*       insert(const std::string& name, const GOptimizerPar& par);
    void                 remove(const int& index);
    void                 remove(const std::string& name);
    void                 reserve(const int& num);
    void                 extend(const GOptimizerPars& pars);
    bool                 contains(const std::string& name) const;
    std::string          print(const GChatter& chatter = NORMAL) const;

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
 * @brief Return pointer to parameter
 *
 * @param[in] index Parameter index [0,...,size()-1].
 * @return Pointer to parameter.
 *
 * Returns a pointer to the parameter with the specified @p index.
 ***************************************************************************/
inline
GOptimizerPar* GOptimizerPars::operator[](const int& index)
{
    return (m_pars[index]);
}


/***********************************************************************//**
 * @brief Return pointer to parameter (const version)
 *
 * @param[in] index Parameter index [0,...,size()-1].
 * @return Pointer to parameter.
 *
 * Returns a pointer to the parameter with the specified @p index.
 ***************************************************************************/
inline
const GOptimizerPar* GOptimizerPars::operator[](const int& index) const
{
    return (m_pars[index]);
}


/***********************************************************************//**
 * @brief Return number of parameters in container
 *
 * @return Number of parameters in container.
 *
 * Returns the number of parameters in the parameter container.
 ***************************************************************************/
inline
int GOptimizerPars::size(void) const
{
    return (m_pars.size());
}


/***********************************************************************//**
 * @brief Signals if there are no parameters in container
 *
 * @return True if container is empty, false otherwise.
 *
 * Signals if the parameters container does not contain any parameter.
 ***************************************************************************/
inline
bool GOptimizerPars::is_empty(void) const
{
    return (m_pars.empty());
}


/***********************************************************************//**
 * @brief Reserves space for parameters in container
 *
 * @param[in] num Number of parameters.
 *
 * Reserves space for @p num parameters in the container.
 ***************************************************************************/
inline
void GOptimizerPars::reserve(const int& num)
{
    m_pars.reserve(num);
    return;
}

#endif /* GOPTIMIZERPARS_HPP */
