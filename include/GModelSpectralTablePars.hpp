/***************************************************************************
 * GModelSpectralTablePars.hpp - Spectral table model par container class  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2019-2020 by Juergen Knoedlseder                         *
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
 * @file GModelSpectralTablePars.hpp
 * @brief Spectral table model parameter container class definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPECTRALTABLEPARS_HPP
#define GMODELSPECTRALTABLEPARS_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GContainer.hpp"
#include "GModelSpectralTablePar.hpp"

/* __ Forward declarations _______________________________________________ */


/***********************************************************************//**
 * @class GModelSpectralTablePars
 *
 * @brief Spectral table model parameter container class
 ***************************************************************************/
class GModelSpectralTablePars : public GContainer {

public:
    // Constructors and destructors
    GModelSpectralTablePars(void);
    GModelSpectralTablePars(const GModelSpectralTablePars& pars);
    virtual ~GModelSpectralTablePars(void);

    // Operators
    GModelSpectralTablePars&      operator=(const GModelSpectralTablePars& pars);
    GModelSpectralTablePar*       operator[](const int& index);
    const GModelSpectralTablePar* operator[](const int& index) const;
    GModelSpectralTablePar*       operator[](const std::string& name);
    const GModelSpectralTablePar* operator[](const std::string& name) const;

    // Implemented pure virtual base class methods
    void                          clear(void);
    GModelSpectralTablePars*      clone(void) const;
    std::string                   classname(void) const;
    GModelSpectralTablePar*       at(const int& index);
    const GModelSpectralTablePar* at(const int& index) const;
    int                           size(void) const;
    bool                          is_empty(void) const;
    GModelSpectralTablePar*       set(const int&                    index,
                                      const GModelSpectralTablePar& par);
    GModelSpectralTablePar*       set(const std::string&            name,
                                      const GModelSpectralTablePar& par);
    GModelSpectralTablePar*       append(const GModelSpectralTablePar& par);
    GModelSpectralTablePar*       insert(const int&                    index,
                                         const GModelSpectralTablePar& par);
    GModelSpectralTablePar*       insert(const std::string&            name,
                                         const GModelSpectralTablePar& par);
    void                          remove(const int& index);
    void                          remove(const std::string& name);
    void                          reserve(const int& num);
    void                          extend(const GModelSpectralTablePars& pars);
    bool                          contains(const std::string& name) const;
    std::string                   print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralTablePars& model);
    void free_members(void);
    int  get_index(const std::string& name) const;

    // Proteced members
    std::vector<GModelSpectralTablePar*> m_pars;  //!< List of parameters
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpectralTablePars").
 ***************************************************************************/
inline
std::string GModelSpectralTablePars::classname(void) const
{
    return ("GModelSpectralTablePars");
}


/***********************************************************************//**
 * @brief Return pointer to table model parameter
 *
 * @param[in] index Table model parameter index [0,...,size()-1].
 *
 * Returns a pointer to the table model parameter with the specified
 * @p index.
 ***************************************************************************/
inline
GModelSpectralTablePar* GModelSpectralTablePars::operator[](const int& index)
{
    return (m_pars[index]);
}


/***********************************************************************//**
 * @brief Return pointer to table model parameter (const version)
 *
 * @param[in] index Table model parameter index [0,...,size()-1].
 *
 * Returns a const pointer to the table model parameter with the specified
 * @p index.
 ***************************************************************************/
inline
const GModelSpectralTablePar* GModelSpectralTablePars::operator[](const int& index) const
{
    return (m_pars[index]);
}


/***********************************************************************//**
 * @brief Return pointer to table model parameter
 *
 * @param[in] name Table model parameter name.
 *
 * Returns a pointer to the table model parameter with the specified @p name.
 ***************************************************************************/
inline
GModelSpectralTablePar* GModelSpectralTablePars::operator[](const std::string& name)
{
    return const_cast<GModelSpectralTablePar*>(static_cast<const GModelSpectralTablePars&>(*this)[name]);
}


/***********************************************************************//**
 * @brief Return number of table model parameters in container
 *
 * @return Number of table model parameters in container.
 *
 * Returns the number of table model parameters in the model container.
 ***************************************************************************/
inline
int GModelSpectralTablePars::size(void) const
{
    return (int)m_pars.size();
}


/***********************************************************************//**
 * @brief Signals if there are no table model parameters in container
 *
 * @return True if container is empty, false otherwise.
 *
 * Signals if the table model parameters container does not contain any table
 * model parameter.
 ***************************************************************************/
inline
bool GModelSpectralTablePars::is_empty(void) const
{
    return (m_pars.empty());
}


/***********************************************************************//**
 * @brief Reserves space for table model parameters in container
 *
 * @param[in] num Number of table model parameters
 *
 * Reserves space for @p num table model parameters in the container.
 ***************************************************************************/
inline
void GModelSpectralTablePars::reserve(const int& num)
{
    m_pars.reserve(num);
    return;
}

#endif /* GMODELSPECTRALTABLEPARS_HPP */
