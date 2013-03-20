/***************************************************************************
 *                    GModelPar.cpp - Model parameter class                *
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
 * @file GModelPar.cpp
 * @brief GModelPar class implementation.
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GModelPar.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_CONSTRUCT    "GModelPar::GModelPar(std::string&, double&, double&)"
#define G_FACTOR_VALUE                     "GModelPar::factor_value(double&)"
#define G_FACTOR_MIN                         "GModelPar::factor_min(double&)"
#define G_FACTOR_MAX                         "GModelPar::factor_max(double&)"
#define G_SCALE                                   "GModelPar::scale(double&)"
#define G_READ                                "GModelPar::read(GXmlElement&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GModelPar::GModelPar(void)
{
    // Initialise members
    init_members();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Parameter constructor
 *
 * @param[in] name Parameter name.
 * @param[in] value Parameter value.
 *
 * Constructs a model parameter from a parameter @p name and a parameter
 * @p value.
 *
 * The parameter is auto-scaled, which for a @p value that differs from zero
 * sets the scale factor to @p value and the @p factor_value to unity. For a
 * @p value of zero, the scale factor will be set to unity and the 
 * @p factor_value will be set to @p value.
 ***************************************************************************/
GModelPar::GModelPar(const std::string& name, const double& value)
{
    // Initialise members
    init_members();

    // Set name attribute
    m_name = name;

    // Set auto-scaled factor and scale attributes
    if (value != 0.0) {
        m_factor_value = 1.0;
        m_scale        = value;
    }
    else {
        m_factor_value = value;
        m_scale        = 1.0;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Parameter constructor
 *
 * @param[in] name Parameter name.
 * @param[in] factor Parameter value factor.
 * @param[in] scale Parameter scaling (non-zero value).
 *
 * @exception GException::invalid_argument
 *            Sacle factor of 0 specified.
 *
 * Constructs a model parameter from a parameter @p name, value @p factor
 * and @p scale factor. The @p scale factor needs to be a non-zero value.
 * If the @p scale factor is zero, an exception is thrown.
 ***************************************************************************/
GModelPar::GModelPar(const std::string& name,
                     const double&      factor,
                     const double&      scale)
{
    // Make sure that scale is not zero
    if (scale == 0.0) {
        std::string msg = "Specified a model scale factor of 0.\n"
                          "Model parameters need a non-zero scale factor.";
        throw GException::invalid_argument(G_CONSTRUCT, msg);
    }

    // Initialise members
    init_members();

    // Set attributes
    m_name         = name;
    m_factor_value = factor;
    m_scale        = scale;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] par Model parameter.
 ***************************************************************************/
GModelPar::GModelPar(const GModelPar& par)
{ 
    // Initialise members
    init_members();

    // Copy members
    copy_members(par);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GModelPar::~GModelPar(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] par Model parameter.
 * @return Model parameter.
 ***************************************************************************/
GModelPar& GModelPar::operator=(const GModelPar& par)
{ 
    // Execute only if object is not identical
    if (this != &par) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(par);

    } // endif: object was not identical
  
    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear model parameter
 *
 * Resets model parameter to a clean initial state.
 ***************************************************************************/
void GModelPar::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone model parameter
 *
 * @return Pointer to deep copy of model parameter.
 ***************************************************************************/
GModelPar* GModelPar::clone(void) const
{
    // Clone model parameter
    return new GModelPar(*this);
}


/***********************************************************************//**
 * @brief Set parameter value
 *
 * @param[in] value Parameter value.
 *
 * Sets the parameter value. The method stores the value factor which is
 * obtained by dividing the @p value by the scale factor.
 *
 * The method calls factor_value() for assigning the value factor, and this
 * method will verify that the value lies within the specified boundaries.
 ***************************************************************************/
void GModelPar::value(const double& value)
{
    // Set value factor. The GModelPar class makes sure that m_scale is
    // never 0, so no test is needed here
    factor_value(value / m_scale);
	
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set parameter error
 *
 * @param[in] error Parameter error.
 *
 * Sets the parameter error. The method stores the error factor which is
 * obtained by dividing the @p error by the scale factor.
 ***************************************************************************/
void GModelPar::error(const double& error)
{
    // Set error factor. The GModelPar class makes sure that m_scale is
    // never 0, so no test is needed here
    factor_error(error / m_scale);
	
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set parameter gradient
 *
 * @param[in] gradient Parameter gradient.
 *
 * Sets the parameter gradient. The method stores the gradient factor which
 * is obtained by dividing @p gradient by the scale factor.
 ***************************************************************************/
void GModelPar::gradient(const double& gradient)
{
    // Set gradient factor. The GModelPar class makes sure that m_scale is
    // never 0, so no test is needed here
    factor_gradient(gradient / m_scale);
	
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set minimum parameter boundary
 *
 * @param[in] min Parameter minimum.
 *
 * Sets the minimum parameter boundary. The method stores the minimum
 * boundary factor which is obtained by dividing @p min by the scale factor.
 ***************************************************************************/
void GModelPar::min(const double& min)
{
    // Set minimum boundary factor. The GModelPar class makes sure that
    // m_scale is never 0, so no test is needed here
    factor_min(min / m_scale);
	
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set maximum parameter boundary
 *
 * @param[in] max Parameter maximum.
 *
 * Sets the maximum parameter boundary. The method stores the maximum
 * boundary factor which is obtained by dividing @p max by the scale factor.
 ***************************************************************************/
void GModelPar::max(const double& max)
{
    // Set maximum boundary factor. The GModelPar class makes sure that
    // m_scale is never 0, so no test is needed here
    factor_max(max / m_scale);
	
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set minimum and maximum parameter boundaries
 *
 * @param[in] min Parameter minimum.
 * @param[in] max Parameter maximum.
 *
 * Sets the minimum and maximum parameter boundaries. The method calls the
 * min() and max() methods to set the boundaries.
 ***************************************************************************/
void GModelPar::range(const double& min, const double& max)
{
    // Set minimum and maximum
    this->min(min);
    this->max(max);
	
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set parameter value factor
 *
 * @param[in] value Value factor.
 *
 * @exception GException::invalid_argument
 *            Parameter @p value outside [min,max] boundaries.
 *
 * Sets the value factor of the model parameter. The method makes sure that
 * none of the boundaries is violated. Otherwise, exceptions are thrown.
 ***************************************************************************/
void GModelPar::factor_value(const double& value)
{
    // If there is a minimum boundary and if value is below this boundary
    // then throw an exception
    if (m_hasmin && value < m_factor_min) {
        std::string msg = "Specified value factor "+str(value)+
                          " is smaller than the minimum boundary "+
                          str(m_factor_min)+".";
        throw GException::invalid_argument(G_FACTOR_VALUE, msg);
    }

    // If there is a maximum boundary and if value is above this boundary
    // then throw an exception
    if (m_hasmax && value > m_factor_max) {
        std::string msg = "Specified value factor "+str(value)+
                          " is larger than the maximum boundary "+
                          str(m_factor_max)+".";
        throw GException::invalid_argument(G_FACTOR_VALUE, msg);
    }

    // Assign value
    m_factor_value = value;
	
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set parameter minimum factor
 *
 * @param[in] min Minimum factor.
 *
 * @exception GException::invalid_argument
 *            Parameter @p min larger than value factor.
 *
 * Sets the minimum boundary factor of the model parameter. The method makes
 * sure that the minimum is not larger than the actual value factor.
 * Otherwise, an exception is thrown.
 ***************************************************************************/
void GModelPar::factor_min(const double& min)
{
    // Check if minimum is larger than value
    if (min > m_factor_value) {
        std::string msg = "Specified minimum factor "+str(min)+
                          " is larger than the value factor "+
                          str(m_factor_value)+".";
        throw GException::invalid_argument(G_FACTOR_VALUE, msg);
    }

    // Assign minimum factor
    m_factor_min = min;
    
    // Flag that minimum was set
    m_hasmin = true;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set parameter maximum factor
 *
 * @param[in] max Maximum factor.
 *
 * @exception GException::invalid_argument
 *            Parameter @p max smaller than value factor.
 *
 * Sets the maximum boundary factor of the model parameter. The method makes
 * sure that the maximum is not smaller than the actual value factor.
 * Otherwise, an exception is thrown.
 ***************************************************************************/
void GModelPar::factor_max(const double& max)
{
    // Check if maximum is smaller than value
    if (max < m_factor_value) {
        std::string msg = "Specified maximum factor "+str(max)+
                          " is smaller than the value factor "+
                          str(m_factor_value)+".";
        throw GException::invalid_argument(G_FACTOR_VALUE, msg);
    }

    // Assign maximum
    m_factor_max = max;
    
    // Flag that maximum was set
    m_hasmax = true;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set parameter minimum and maximum factors
 *
 * @param[in] min Minimum factor.
 * @param[in] max Maximum factor.
 *
 * Sets the minimum and maximum boundary factors. The method calls the
 * factor_min() and factor_max() methods which perform validity checking
 * of the arguments.
 ***************************************************************************/
void GModelPar::factor_range(const double& min, const double& max)
{
    // Set minimum and maximum
    factor_min(min);
    factor_max(max);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set scale factor
 *
 * @param[in] scale Scale factor.
 *
 * @exception GException::invalid_argument
 *            Sacle factor of 0 specified.
 *
 * Sets the scale factor of the model parameter. All parameter attributes are
 * rescaled accordingly.
 *
 * Special care is taken for value boundaries, that need swaping if the scale
 * factor changes its sign. In case of sign change the following logic applies:
 * - if minimum and maximum boundaries are set they are swapped
 * - if only a minimum boundary exists it will be replaced by a maximum boundary
 * - if only a maximum boundary exists it will be replaced by a minimum boundary
 *
 * An exception is thrown if a scale factor of 0 is specified.
 ***************************************************************************/
void GModelPar::scale(const double& scale)
{
    // Make sure that scale is not zero
    if (scale == 0.0) {
        std::string msg = "Specified scale factor of 0.\n"
                          "Model parameters need a non-zero scale factor.";
        throw GException::invalid_argument(G_SCALE, msg);
    }

    // Set rescaling
    double rescale = m_scale/scale;

    // Set new scale factor
    m_scale = scale;

    // Set values, error, gradient, min and max
    m_factor_value    *= rescale;
    m_factor_error    *= rescale;
    m_factor_gradient *= rescale;
    if (m_hasmin) {
        m_factor_min *= rescale;
    }
    if (m_hasmax) {
        m_factor_max *= rescale;
    }

    // Takes care of boundaries in case of sign change
    if (rescale < 0.0) {
        if (m_hasmin && m_hasmax) {
            double swap  = m_factor_min;
            m_factor_min = m_factor_max;
            m_factor_max = swap;
        }
        else if (m_hasmin) {
            m_factor_max = m_factor_min;
            m_hasmax     = true;
            m_factor_min = 0.0;
            m_hasmin     = false;
        }
        else if (m_hasmax) {
            m_factor_min = m_factor_max;
            m_hasmin     = true;
            m_factor_max = 0.0;
            m_hasmax     = false;
        }
    }
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Autoscale parameter
 *
 * Sets the value factor to unity and the scale factor to the real value
 * of the parameter. The method will also adjust the error factor and
 * gradient factor, as well as the minimum and maximum factors if they
 * exist.
 *
 * Special care is taken for value boundaries, that need swaping if the scale
 * factor changes its sign. In case of sign change the following logic applies:
 * - if minimum and maximum boundaries are set they are swapped
 * - if only a minimum boundary exists it will be replaced by a maximum boundary
 * - if only a maximum boundary exists it will be replaced by a minimum boundary
 *
 * The method does nothing if the actual value factor is zero.
 ***************************************************************************/
void GModelPar::autoscale(void)
{
    // Continue only if the value factor is non-zero
    if (m_factor_value != 0.0) {

        // Set the scale factor
        m_scale *= m_factor_value;

        // Get inverse scaling factor
        double invscale = 1.0 / m_factor_value;

        // Set values, error, gradient, min and max
        m_factor_value    *= invscale;
        m_factor_error    *= invscale;
        m_factor_gradient *= invscale;
        if (m_hasmin) {
            m_factor_min *= invscale;
        }
        if (m_hasmax) {
            m_factor_max *= invscale;
        }

        // Takes care of boundaries in case of sign change
        if (invscale < 0.0) {
            if (m_hasmin && m_hasmax) {
                double swap  = m_factor_min;
                m_factor_min = m_factor_max;
                m_factor_max = swap;
            }
            else if (m_hasmin) {
                m_factor_max = m_factor_min;
                m_hasmax     = true;
                m_factor_min = 0.0;
                m_hasmin     = false;
            }
            else if (m_hasmax) {
                m_factor_min = m_factor_max;
                m_hasmin     = true;
                m_factor_max = 0.0;
                m_hasmax     = false;
            }
        }

    } // endif: value was non-zero

    // Return
    return;
}


/***********************************************************************//**
 * @brief Extract parameter attributes from XML element
 *
 * @param[in] xml XML element
 *
 * @exception GException::invalid_value
 *            Invalid combination of parameter attributes encountered.
 *
 * Extracts the parameter attributes from an XML element of the form
 *
 *     <parameter name=".." value=".." error=".." scale=".." min=".." max="..' free="..">
 *
 * Each of the attributes are optional, with the following scheme for
 * assigning default values in case that the attribute was not found:
 * - @p value sets @p m_factor_value (defaults to 0.0)
 * - @p error sets @p m_factor_error (defaults to 0.0)
 * - @p scale sets @p m_scale (defaults to 1.0)
 * - @p min sets @p m_factor_min (will remove_min() if not found)
 * - @p max sets @p m_factor_max (will remove_max() if not found)
 * - @p free sets @p m_free (papameter will be fixed if not found)
 ***************************************************************************/
void GModelPar::read(const GXmlElement& xml)
{
    // Get value
    std::string arg = xml.attribute("value");
    if (arg != "") {
        m_factor_value = todouble(arg);
    }
    else {
        m_factor_value = 0.0;
    }

    // Get error
    arg = xml.attribute("error");
    if (arg != "") {
        m_factor_error = todouble(arg);
    }
    else {
        m_factor_error = 0.0;
    }

    // Get scale factor
    arg = xml.attribute("scale");
    if (arg != "") {
        m_scale = todouble(arg);
    }
    else {
        m_scale = 1.0;
    }

    // Get min
    arg = xml.attribute("min");
    if (arg != "") {
        m_factor_min = todouble(arg);
        m_hasmin     = true;
    }
    else {
        remove_min();
    }

    // Get max
    arg = xml.attribute("max");
    if (arg != "") {
        m_factor_max = todouble(arg);
        m_hasmax     = true;
    }
    else {
        remove_max();
    }

    // Get free
    if (xml.attribute("free") == "1" || 
        tolower(xml.attribute("free")) == "true") {
        free();
    }
    else {
        fix();
    }

    // If there is a minimum and maximum, make sure that the maximum is
    // not smaller than the minimum
    if (m_hasmin && m_hasmax) {
        if (m_factor_min > m_factor_max) {
            std::string msg = "The model parameter "+m_name+
                              " in the XML document has a minimum boundary "+
                              str(m_factor_min)+
                              " that is larger than the maximum boundary "+
                              str(m_factor_max)+".\n"+xml.print();
            throw GException::invalid_value(G_READ, msg);
        }
    }

    // If there is a minimum, make sure that the value is not below it
    if (m_hasmin && m_factor_value < m_factor_min) {
        std::string msg = "The model parameter "+m_name+
                          " in the XML document has a value "+
                            str(m_factor_value)+
                            " that is smaller than the minimum boundary "+
                            str(m_factor_min)+".\n"+xml.print();
        throw GException::invalid_value(G_READ, msg);
    }

    // If there is a maximum, make sure that the value is not above it
    if (m_hasmax && m_factor_value > m_factor_max) {
        std::string msg = "The model parameter "+m_name+
                          " in the XML document has a value "+
                            str(m_factor_value)+
                            " that is larger than the maximum boundary "+
                            str(m_factor_max)+".\n"+xml.print();
        throw GException::invalid_value(G_READ, msg);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set or update parameter attributes in XML element
 *
 * @param[in] xml XML element.
 *
 * Sets or updates the parameter attributes in an XML element of the form
 *
 *     <parameter name=".." value=".." error=".." scale=".." min=".." max="..' free="..">
 *
 * The following attributes will be set:
 * - @p value
 * - @p error (only in case that the parameter is free)
 * - @p scale
 * - @p min (only in case that a minimum exists)
 * - @p max (only in case that a maximum exists)
 * - @p free
 ***************************************************************************/
void GModelPar::write(GXmlElement& xml) const
{
    // Set value
    xml.attribute("value", str(m_factor_value));

    // Set error (only if parameter is free)
    if (isfree()) {
        xml.attribute("error", str(m_factor_error));
    }

    // Set scale
    xml.attribute("scale", str(m_scale));

    // Set minimum
    if (hasmin()) {
        xml.attribute("min", str(m_factor_min));
    }

    // Set maximum
    if (hasmax()) {
        xml.attribute("max", str(m_factor_max));
    }

    // Set free/fix flag
    if (isfree()) {
        xml.attribute("free", "1");
    }
    else {
        xml.attribute("free", "0");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print parameter information
 *
 * @return String with model parameter information.
 ***************************************************************************/
std::string GModelPar::print(void) const
{
    // Initialise result string
    std::string result;

    // Append parameter name
    result.append(parformat(" "+name()));

    // Append value
    result.append(str(value()));

    // For free parameters, append statistical uncertainty
    if (m_free) {
        result.append(" +/- "+str(std::abs(error())));
    }

    // Append parameter limites if they exist
    if (m_hasmin && m_hasmax) {
        result.append(" ["+str(min()) + ","+str(max())+"]");
    }
    else if (m_hasmin) {
        result.append(" ["+str(min()) + ",infty[");
    }
    else if (m_hasmax) {
        result.append(" ]-infty,"+str(max())+"]");
    }

    // Append parameter unit
    result.append(" "+m_unit);

    // Signal if parameter was free or fixed
    if (m_free) {
        result.append(" (free");
    }
    else {
        result.append(" (fixed");
    }

    // Append parameter scale
    result.append(",scale="+str(m_scale));

    // Signal if parameter has analytic gradient
    if (m_hasgrad) {
        result.append(",gradient)");
    }
    else {
        result.append(")");
    }

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModelPar::init_members(void)
{
    // Initialise members
    m_name.clear();
    m_unit.clear();
    m_factor_value    = 0.0;
    m_factor_error    = 0.0;
    m_factor_gradient = 0.0;
    m_factor_min      = 0.0;
    m_factor_max      = 0.0;
    m_scale           = 1.0;
    m_free            = true;
    m_hasmin          = false;
    m_hasmax          = false;
    m_hasgrad         = false;
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] par Model parameter.
 ***************************************************************************/
void GModelPar::copy_members(const GModelPar& par)
{
    // Copy members
    m_name            = par.m_name;
    m_unit            = par.m_unit;
    m_factor_value    = par.m_factor_value;
    m_factor_error    = par.m_factor_error;
    m_factor_gradient = par.m_factor_gradient;
    m_factor_min      = par.m_factor_min;
    m_factor_max      = par.m_factor_max;
    m_scale           = par.m_scale;
    m_free            = par.m_free;
    m_hasmin          = par.m_hasmin;
    m_hasmax          = par.m_hasmax;
    m_hasgrad         = par.m_hasgrad;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelPar::free_members(void)
{  
    // Return
    return;
}
