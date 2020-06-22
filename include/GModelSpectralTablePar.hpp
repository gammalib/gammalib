/***************************************************************************
 *    GModelSpectralTablePar.hpp - Spectral table model parameter class    *
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
 * @file GModelSpectralTablePar.hpp
 * @brief Spectral table model parameter class definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPECTRALTABLEPAR_HPP
#define GMODELSPECTRALTABLEPAR_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GModelPar.hpp"
#include "GNodeArray.hpp"

/* __ Forward declarations _______________________________________________ */


/***********************************************************************//**
 * @class GModelSpectralTablePar
 *
 * @brief Spectral table model parameter class
 ***************************************************************************/
class GModelSpectralTablePar : public GBase {

public:
    // Constructors and destructors
    GModelSpectralTablePar(void);
    GModelSpectralTablePar(const GModelPar&           par,
                           const std::vector<double>& values);
    GModelSpectralTablePar(const GModelSpectralTablePar& par);
    virtual ~GModelSpectralTablePar(void);

    // Operators
    GModelSpectralTablePar& operator=(const GModelSpectralTablePar& par);

    // Methods
    void                    clear(void);
    GModelSpectralTablePar* clone(void) const;
    std::string             classname(void) const;
    int                     size(void) const;
    bool                    is_empty(void) const;
    GModelPar&              par(void);
    const GModelPar&        par(void) const;
    const GNodeArray&       values(void) const;
    const int&              method(void) const;
    void                    method(const int& method);
    std::string             print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralTablePar& par);
    void free_members(void);

    // Protected members
    GModelPar  m_par;    //!< Model parameter
    GNodeArray m_values; //!< Parameter values
    int        m_method; //!< Interpolation method (0: linear, 1: logarithmic)
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpectralTablePar").
 ***************************************************************************/
inline
std::string GModelSpectralTablePar::classname(void) const
{
    return ("GModelSpectralTablePar");
}


/***********************************************************************//**
 * @brief Return number of table model parameter values
 *
 * @return Number of table model parameter values.
 *
 * Returns the number of table model parameter values.
 ***************************************************************************/
inline
int GModelSpectralTablePar::size(void) const
{
    return (m_values.size());
}


/***********************************************************************//**
 * @brief Signals if there are no table model parameter values
 *
 * @return True if there are no table model parameter values, false otherwise.
 *
 * Signals if there are no table model parameter values.
 ***************************************************************************/
inline
bool GModelSpectralTablePar::is_empty(void) const
{
    return (m_values.is_empty());
}


/***********************************************************************//**
 * @brief Return reference to table model parameter
 *
 * @return Reference to table model parameter.
 ***************************************************************************/
inline
GModelPar& GModelSpectralTablePar::par(void)
{
    return m_par;
}


/***********************************************************************//**
 * @brief Return reference to table model parameter (const version)
 *
 * @return Reference to table model parameter.
 ***************************************************************************/
inline
const GModelPar& GModelSpectralTablePar::par(void) const
{
    return m_par;
}


/***********************************************************************//**
 * @brief Return reference to table model parameter values as node array
 *
 * @return Reference to table model parameter values as node array.
 ***************************************************************************/
inline
const GNodeArray& GModelSpectralTablePar::values(void) const
{
    return m_values;
}


/***********************************************************************//**
 * @brief Return reference to table model parameter interpolation method
 *
 * @return Interpolation method (0: linear, 1: logarithmic).
 ***************************************************************************/
inline
const int& GModelSpectralTablePar::method(void) const
{
    return m_method;
}


/***********************************************************************//**
 * @brief Set table model parameter interpolation method
 *
 * @param[in] method Interpolation method (0: linear, 1: logarithmic).
 ***************************************************************************/
inline
void GModelSpectralTablePar::method(const int& method)
{
    m_method = method;
    return;
}

#endif /* GMODELSPECTRALTABLEPAR_HPP */
