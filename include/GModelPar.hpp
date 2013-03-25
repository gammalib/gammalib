/***************************************************************************
 *                   GModelPar.hpp - Model parameter class                 *
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
 * @file GModelPar.hpp
 * @brief Model parameter class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELPAR_HPP
#define GMODELPAR_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
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
 * attributes is tested using the hasmin() and hasmax() methods,
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
 * methods and it is retrieved using the isfree() and isfixed() methods.
 * The attribute @p grad is set and retrieved using the has hasgrad()
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
class GModelPar : public GBase {

public:
    // Constructors and destructors
    GModelPar(void);
    explicit GModelPar(const std::string& name, const double& value);
    explicit GModelPar(const std::string& name, const double& factor, const double& scale);
    GModelPar(const GModelPar& par);
    virtual ~GModelPar(void);

    // Operators
    GModelPar& operator=(const GModelPar& par);

    // Attribute methods 
    double value(void) const;
    double error(void) const;
    double gradient(void) const;
    double min(void) const;
    double max(void) const;
    void   value(const double& value);
    void   error(const double& error);
    void   gradient(const double& gradient);
    void   min(const double& min);
    void   max(const double& max);
    void   range(const double& min, const double& max);

    // Factorization methods
    const double& factor_value(void) const;
    const double& factor_error(void) const;
    const double& factor_gradient(void) const;
    const double& factor_min(void) const;
    const double& factor_max(void) const;
    const double& scale(void) const;
    void          factor_value(const double& value);
    void          factor_error(const double& error);
    void          factor_gradient(const double& gradient);
    void          factor_min(const double& min);
    void          factor_max(const double& max);
    void          factor_range(const double& min, const double& max);
    void          scale(const double& scale);

    // Boundary methods
    bool hasmin(void) const;
    bool hasmax(void) const;
    bool hasrange(void) const;
    void remove_min(void);
    void remove_max(void);
    void remove_range(void);

    // Property methods
    bool isfree(void) const;
    bool isfixed(void) const;
    bool hasgrad(void) const;
    void free(void);
    void fix(void);
    void hasgrad(const bool& grad);

    // Other methods
    void               clear(void);
    GModelPar*         clone(void) const;
    const std::string& name(void) const;
    const std::string& unit(void) const;
    void               name(const std::string& name);
    void               unit(const std::string& unit);
    void               autoscale(void);
    void               read(const GXmlElement& xml);
    void               write(GXmlElement& xml) const;
    std::string        print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelPar& par);
    void free_members(void);

    // Proteced data members
    std::string m_name;            //!< Parameter name
    std::string m_unit;            //!< Parameter unit
    double      m_factor_value;    //!< Parameter value factor
    double      m_factor_error;    //!< Uncertainty in parameter value factor
    double      m_factor_gradient; //!< Model gradient factor
    double      m_factor_min;      //!< Parameter minimum factor
    double      m_factor_max;      //!< Parameter maximum factor
    double      m_scale;           //!< Parameter scaling (true = factor * scale)
    bool        m_free;            //!< Parameter is free
    bool        m_hasmin;          //!< Parameter has minimum boundary
    bool        m_hasmax;          //!< Parameter has maximum boundary
    bool        m_hasgrad;         //!< Parameter has analytic gradient
};


/***********************************************************************//**
 * @brief Return parameter value
 *
 * @return Parameter value.
 *
 * Returns the parameter value. The parameter value is computed by
 * multiplying the value factor by the scale factor.
 ***************************************************************************/
inline
double GModelPar::value(void) const
{
    return (m_factor_value * m_scale);
}


/***********************************************************************//**
 * @brief Return parameter error
 *
 * @return Parameter error.
 *
 * Returns the parameter error. The parameter error is computed by
 * multiplying the error factor by the scale factor.
 ***************************************************************************/
inline
double GModelPar::error(void) const
{
    return (m_factor_error * m_scale);
}


/***********************************************************************//**
 * @brief Return parameter gradient
 *
 * @return Parameter gradient.
 *
 * Returns the parameter gradient. The parameter gradient is computed by
 * multiplying the gradient factor by the scale factor.
 ***************************************************************************/
inline
double GModelPar::gradient(void) const
{
    return (m_factor_gradient * m_scale);
}


/***********************************************************************//**
 * @brief Return parameter minimum boundary
 *
 * @return Parameter minimum boundary.
 *
 * Returns the parameter minimum boundary. The parameter minimum boundary is
 * computed by multiplying the minimum boundary factor by the scale factor.
 ***************************************************************************/
inline
double GModelPar::min(void) const
{
    return (m_factor_min * m_scale);
}


/***********************************************************************//**
 * @brief Return parameter maximum boundary
 *
 * @return Parameter maximum boundary.
 *
 * Returns the parameter maximum boundary. The parameter maximum boundary is
 * computed by multiplying the maximum boundary factor by the scale factor.
 ***************************************************************************/
inline
double GModelPar::max(void) const
{
    return (m_factor_max * m_scale);
}


/***********************************************************************//**
 * @brief Return parameter value factor
 *
 * @return Parameter value factor.
 *
 * Returns the parameter value factor.
 ***************************************************************************/
inline
const double& GModelPar::factor_value(void) const
{
    return (m_factor_value);
}


/***********************************************************************//**
 * @brief Return parameter error factor
 *
 * @return Parameter error factor.
 *
 * Returns the parameter error factor.
 ***************************************************************************/
inline
const double& GModelPar::factor_error(void) const
{
    return (m_factor_error);
}


/***********************************************************************//**
 * @brief Return parameter gradient factor
 *
 * @return Parameter gradient factor.
 *
 * Returns the parameter gradient factor.
 ***************************************************************************/
inline
const double& GModelPar::factor_gradient(void) const
{
    return (m_factor_gradient);
}


/***********************************************************************//**
 * @brief Return parameter minimum boundary factor
 *
 * @return Parameter minimum boundary factor.
 *
 * Returns the parameter minimum boundary factor.
 ***************************************************************************/
inline
const double& GModelPar::factor_min(void) const
{
    return (m_factor_min);
}


/***********************************************************************//**
 * @brief Return parameter maximum boundary factor
 *
 * @return Parameter maximum boundary factor.
 *
 * Returns the parameter maximum boundary factor.
 ***************************************************************************/
inline
const double& GModelPar::factor_max(void) const
{
    return (m_factor_max);
}


/***********************************************************************//**
 * @brief Return parameter scale
 *
 * @return Parameter scale factor.
 *
 * Returns the parameter scale factor.
 ***************************************************************************/
inline
const double& GModelPar::scale(void) const
{
    return (m_scale);
}


/***********************************************************************//**
 * @brief Set parameter error factor
 *
 * @param[in] error Parameter error factor.
 *
 * Sets the parameter error factor.
 ***************************************************************************/
inline
void GModelPar::factor_error(const double& error)
{
    m_factor_error = error;
    return;
}


/***********************************************************************//**
 * @brief Set parameter gradient factor
 *
 * @param[in] gradient Parameter gradient factor.
 *
 * Sets the parameter gradient factor.
 ***************************************************************************/
inline
void GModelPar::factor_gradient(const double& gradient)
{
    m_factor_gradient = gradient;
    return;
}


/***********************************************************************//**
 * @brief Signal if parameter has minimum boundary
 *
 * @return True if parameter has minimum boundary, false otherwise.
 *
 * Signals if the parameter has a minimum boundary.
 ***************************************************************************/
inline
bool GModelPar::hasmin(void) const
{
    return m_hasmin;
}


/***********************************************************************//**
 * @brief Signal if parameter has maximum boundary
 *
 * @return True if parameter has maximum boundary, false otherwise.
 *
 * Signals if the parameter has a maximum boundary.
 ***************************************************************************/
inline
bool GModelPar::hasmax(void) const
{
    return m_hasmax;
}


/***********************************************************************//**
 * @brief Signal if parameter has minimum and maximum boundaries
 *
 * @return True if parameter has minimum and maximum boundaries, false
 *         otherwise.
 *
 * Signals if the parameter has a minimum and a maximum boundary.
 ***************************************************************************/
inline
bool GModelPar::hasrange(void) const
{
    return (m_hasmin && m_hasmax);
}


/***********************************************************************//**
 * @brief Removes minimum boundary
 *
 * Removes minimum boundary from the model parameter.
 ***************************************************************************/
inline
void GModelPar::remove_min(void)
{
    m_hasmin = false;
    return;
}


/***********************************************************************//**
 * @brief Removes maximum boundary
 *
 * Removes maximum boundary from the model parameter.
 ***************************************************************************/
inline
void GModelPar::remove_max(void)
{
    m_hasmax = false;
    return;
}


/***********************************************************************//**
 * @brief Removes minimum and maximum boundary
 *
 * Removes minimum and maximum boundary from the model parameter.
 ***************************************************************************/
inline
void GModelPar::remove_range(void)
{
    m_hasmin = false;
    m_hasmax = false;
    return;
}


/***********************************************************************//**
 * @brief Signal if parameter is free
 *
 * @return True if parameter is free, false otherwise.
 *
 * Signals if the parameter is free, i.e. that it shall be fitted in a
 * parameter optimization process.
 ***************************************************************************/
inline
bool GModelPar::isfree(void) const
{
    return m_free;
}


/***********************************************************************//**
 * @brief Signal if parameter is fixed
 *
 * @return True if parameter is fixed, false otherwise.
 *
 * Signals if the parameter is fixed, i.e. that it shall NOT be fitted in a
 * parameter optimization process.
 ***************************************************************************/
inline
bool GModelPar::isfixed(void) const
{
    return (!m_free);
}


/***********************************************************************//**
 * @brief Signal if parameter gradient is computed analytically
 *
 * @return True if parameter is gradient is computed analytically, false
 *         otherwise.
 *
 * Signals if the parameter gradient is computed analytically. This property
 * is used in the model optimization process to identify parameters for which
 * gradients need to be computed numerically.
 ***************************************************************************/
inline
bool GModelPar::hasgrad(void) const
{
    return m_hasgrad;
}


/***********************************************************************//**
 * @brief Free a parameter
 *
 * Frees a parameter for model fitting. The parameter shall be fitted in a
 * parameter optimization process.
 ***************************************************************************/
inline
void GModelPar::free(void)
{
    m_free = true;
    return;
}


/***********************************************************************//**
 * @brief Fix a parameter
 *
 * Fixes a parameter for model fitting. The parameter shall NOT be altered
 * in a parameter optimization process.
 ***************************************************************************/
inline
void GModelPar::fix(void)
{
    m_free = false;
    return;
}


/***********************************************************************//**
 * @brief Set gradient property
 *
 * @param[in] grad Gradient flag
 *
 * Sets the gradient property of the model parameter. If @p grad is set to
 * true, the parameter gradient will be computed analytically by the method
 * that evaluates the model. If false, the gradients needs to be computed
 * numerically.
 ***************************************************************************/
inline
void GModelPar::hasgrad(const bool& grad)
{
    m_hasgrad = grad;
    return;
}


/***********************************************************************//**
 * @brief Return parameter name
 *
 * @return Parameter name.
 *
 * Returns the parameter name.
 ***************************************************************************/
inline
const std::string& GModelPar::name(void) const
{
    return m_name;
}


/***********************************************************************//**
 * @brief Return parameter unit
 *
 * @return Parameter unit.
 *
 * Returns the parameter unit.
 ***************************************************************************/
inline
const std::string& GModelPar::unit(void) const
{
    return m_unit;
}


/***********************************************************************//**
 * @brief Set parameter name
 *
 * @param[in] name Parameter name.
 *
 * Sets the parameter name.
 ***************************************************************************/
inline
void GModelPar::name(const std::string& name)
{
    m_name = name;
    return;
}


/***********************************************************************//**
 * @brief Set parameter unit
 *
 * @param[in] unit Parameter unit.
 *
 * Sets the parameter unit.
 ***************************************************************************/
inline
void GModelPar::unit(const std::string& unit)
{
    m_unit = unit;
    return;
}

#endif /* GMODELPAR_HPP */
