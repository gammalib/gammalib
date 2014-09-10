/***************************************************************************
 *                   GModelPar.hpp - Model parameter class                 *
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
 * @file GModelPar.hpp
 * @brief Model parameter class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELPAR_HPP
#define GMODELPAR_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GOptimizerPar.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelPar
 *
 * @brief Model parameter class
 *
 * This class implements a model parameter. A model parameter is a numerical
 * value that is used to describe a model. A model parameter has the
 * following attributes:
 * - @p value gives the numerical value of the parameter
 * - @p error gives the statistical uncertainty in the parameter value
 * - @p gradient gives the gradient of a the model with respect to the
 *      parameter
 * - @p min gives the minimum value that the parameter value can take
 * - @p max gives the maximum value that the parameter value can take
 *
 * The parameter attributes are set and retrieved using the value(),
 * error(), gradient(), min() and max() methods, respectively. Furthermore,
 * the range() method is used to simultaneously set the minimum and maximum
 * value of a parameter.
 *
 * The minimum and maximum values are optional, and existence of these
 * attributes is tested using the has_min() and has_max() methods,
 * respectively. The minimum value, maximum value are removed using the
 * remove_min() and remove_max() methods. Simultaneous removal of minimum
 * and maximum values is done using the remove_range() method.
 *
 * Each parameter has furthermore the following properties:
 * - @p free specifies whether the parameter should be fitted
 * - @p grad specifies whether the parameter gradient is computed
 *      analytically (true) or numerically (false)
 *
 * The parameter property @p free is set using the free() and fix()
 * methods and it is retrieved using the is_free() and is_fixed() methods.
 * The attribute @p grad is set and retrieved using the has has_grad()
 * methods.
 *
 * Each model parameter is factorized into a @p factor and a @p scale
 * term. The GModelPar class stores the factors and the scale factor has
 * data members, and the true values are computed using the following
 * relations:
 *
 *     value    = m_factor_value    * m_scale
 *     error    = m_factor_error    * m_scale
 *     gradient = m_factor_gradient * m_scale
 *     min      = m_factor_min      * m_scale
 *     max      = m_factor_max      * m_scale
 *
 * The @p factor and @p scale terms can be set and retrieved using the
 * factor_value(), factor_error(), factor_gradient(), factor_min(),
 * factor_max() and scale() methods.
 ***************************************************************************/
class GModelPar : public GOptimizerPar {

public:
    // Constructors and destructors
    GModelPar(void);
    GModelPar(const std::string& name, const double& value);
    GModelPar(const std::string& name, const double& factor,
              const double& scale);
    GModelPar(const GModelPar& par);
    virtual ~GModelPar(void);

    // Operators
    GModelPar& operator=(const GModelPar& par);

    // Methods
    GModelPar*  clone(void) const;
    std::string classname(void) const;
    void        read(const GXmlElement& xml);
    void        write(GXmlElement& xml) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelPar& par);
    void free_members(void);
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelPar").
 ***************************************************************************/
inline
std::string GModelPar::classname(void) const
{
    return ("GModelPar");
}

#endif /* GMODELPAR_HPP */
