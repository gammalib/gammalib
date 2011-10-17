/***************************************************************************
 *                  GModelPar.hpp  -  Model parameter class                *
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
 * @file GModelPar.hpp
 * @brief Model parameter class interface definition
 * @author J. Knodlseder
 */

#ifndef GMODELPAR_HPP
#define GMODELPAR_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GLog.hpp"
#include "GException.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelPar
 *
 * @brief Model parameter class
 *
 * This method implements a model parameter. The model parameter is
 * factorised into a "value" and a "scale" factor, the real parameter value
 * being the product of both. Parameter optimization is only done on the
 * "value" part of the parameter. This allows scaling of the parameter to
 * have a "value" around unity (which makes the optimizer routines better
 * behave). Methods that work on the real value (i.e. the product of
 * "value" and "scale") are prefixed with "real_".
 *
 * The parameter "value" can optionally by bounded by a minimum and/or
 * maximum value. Note that the boundaries apply to "value".
 *
 * Each parameter has a name, and holds optionally information about its
 * gradient (if "hasgrad" is true).
 ***************************************************************************/
class GModelPar {

    // I/O friends
    friend std::ostream& operator<<(std::ostream& os, const GModelPar& par);
    friend GLog&         operator<<(GLog& log,        const GModelPar& par);

public:
    // Constructors and destructors
    GModelPar(void);
    GModelPar(const GModelPar& par);
    virtual ~GModelPar(void);

    // Operators
    GModelPar& operator=(const GModelPar& par);

    // Methods
    void        clear(void);
    GModelPar*  clone(void) const;
    std::string name(void) const { return m_name; }
    std::string unit(void) const { return m_unit; }
    double      real_value(void) const { return m_value*m_scale; }
    double      real_error(void) const { return m_error*m_scale; }
    double      real_gradient(void) const { return m_gradient*m_scale; }
    double      real_min(void) const { return m_min*m_scale; }
    double      real_max(void) const { return m_max*m_scale; }
    double      value(void) const { return m_value; }
    double      error(void) const { return m_error; }
    double      gradient(void) const { return m_gradient; }
    double      min(void) const { return m_min; }
    double      max(void) const { return m_max; }
    double      scale(void) const { return m_scale; }
    bool        isfree(void) const { return m_free; }
    bool        hasmin(void) const { return m_hasmin; }
    bool        hasmax(void) const { return m_hasmax; }
    bool        hasgrad(void) const { return m_hasgrad; }
    void        name(const std::string& name) { m_name=name; }
    void        unit(const std::string& unit) { m_unit=unit; }
    void        real_value(const double& value);
    void        real_error(const double& error);
    void        value(const double& value);
    void        error(const double& error)  { m_error=error; }
    void        gradient(const double& gradient) { m_gradient=gradient; }
    void        min(const double& min);
    void        max(const double& max);
    void        scale(const double& scale) { m_scale=scale; }
    void        range(const double& min, const double& max);
    void        remove_min(void) { m_hasmin=false; }
    void        remove_max(void) { m_hasmax=false; }
    void        remove_range(void) { m_hasmin=false; m_hasmax=false; }
    void        free(void) { m_free=true; }
    void        fix(void) { m_free=false; }
    void        hasgrad(const bool& grad) { m_hasgrad=grad; }
    void        read(const GXmlElement& xml);
    void        write(GXmlElement& xml) const;
    std::string print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelPar& par);
    void free_members(void);

    // Proteced data members
    std::string  m_name;         //!< Parameter name
    std::string  m_unit;         //!< Parameter unit
    double       m_value;        //!< Parameter value
    double       m_error;        //!< Uncertainty in parameter value
    double       m_gradient;     //!< Model gradient
    double       m_min;          //!< Parameter minimum
    double       m_max;          //!< Parameter maximum
    double       m_scale;        //!< Parameter scale (real = m_value * m_scale)
    bool         m_free;         //!< Parameter is free
    bool         m_hasmin;       //!< Parameter has minimum boundary
    bool         m_hasmax;       //!< Parameter has maximum boundary
    bool         m_hasgrad;      //!< Parameter has analytic gradient
};

#endif /* GMODELPAR_HPP */
