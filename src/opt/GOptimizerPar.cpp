/***************************************************************************
 *                GOptimizerPar.cpp - Optimizer parameter class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2020 by Juergen Knoedlseder                         *
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
 * @file GOptimizerPar.cpp
 * @brief Optimizer parameter class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GOptimizerPar.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_CONSTRUCT     "GOptimizerPar::GOptimizerPar(std::string&, double&,"\
                                                                  " double&)"
#define G_MIN                                   "GOptimizerPar::min(double&)"
#define G_MAX                                   "GOptimizerPar::max(double&)"
#define G_FACTOR_VALUE                 "GOptimizerPar::factor_value(double&)"
#define G_FACTOR_MIN                     "GOptimizerPar::factor_min(double&)"
#define G_FACTOR_MAX                     "GOptimizerPar::factor_max(double&)"
#define G_SCALE                               "GOptimizerPar::scale(double&)"

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
GOptimizerPar::GOptimizerPar(void)
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
 * Constructs a parameter from a parameter @p name and a parameter @p value.
 *
 * The parameter is auto-scaled, which for a @p value that differs from zero
 * sets the scale factor to @p value and the @p factor_value to unity. For a
 * @p value of zero, the scale factor will be set to unity and the 
 * @p factor_value will be set to @p value.
 ***************************************************************************/
GOptimizerPar::GOptimizerPar(const std::string& name, const double& value)
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
 * Constructs a parameter from a parameter @p name, value @p factor
 * and @p scale factor. The @p scale factor needs to be a non-zero value.
 * If the @p scale factor is zero, an exception is thrown.
 ***************************************************************************/
GOptimizerPar::GOptimizerPar(const std::string& name,
                             const double&      factor,
                             const double&      scale)
{
    // Make sure that scale is not zero
    if (scale == 0.0) {
        std::string msg = "Specified a parameter scale factor of 0.\n"
                          "Parameters need a non-zero scale factor.";
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
 * @param[in] par Function parameter.
 ***************************************************************************/
GOptimizerPar::GOptimizerPar(const GOptimizerPar& par)
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
GOptimizerPar::~GOptimizerPar(void)
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
 * @param[in] par Function parameter.
 * @return Function parameter.
 ***************************************************************************/
GOptimizerPar& GOptimizerPar::operator=(const GOptimizerPar& par)
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
 * @brief Clear parameter
 *
 * Resets parameter to a clean initial state.
 ***************************************************************************/
void GOptimizerPar::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone parameter
 *
 * @return Pointer to deep copy of parameter.
 ***************************************************************************/
GOptimizerPar* GOptimizerPar::clone(void) const
{
    // Clone parameter
    return new GOptimizerPar(*this);
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
void GOptimizerPar::value(const double& value)
{
    // Set value factor. The GOptimizerPar class makes sure that m_scale is
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
void GOptimizerPar::error(const double& error)
{
    // Set error factor. The GOptimizerPar class makes sure that m_scale is
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
 * is obtained by multiplying @p gradient by the scale factor.
 ***************************************************************************/
void GOptimizerPar::gradient(const double& gradient)
{
    // Set gradient factor.
    factor_gradient(gradient * m_scale);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set minimum parameter boundary
 *
 * @param[in] min Minimum parameter boundary.
 *
 * @exception GException::invalid_argument
 *            Minimum parameter boundary is larger than the parameter value
 *
 * Sets the minimum parameter boundary.
 ***************************************************************************/
void GOptimizerPar::min(const double& min)
{
    // Check if minimum is larger than value
    if (min > value()) {
        std::string msg = "Specified minimum parameter value "+
                          gammalib::str(min)+" is larger than the parameter "
                          "value "+gammalib::str(value())+".";
        throw GException::invalid_argument(G_MIN, msg);
    }

    // If scale is negative then set maximum factor boundary
    if (m_scale < 0) {
        m_factor_max     = min / m_scale;
        m_has_factor_max = true;
    }

    // ... otherwise set minimum factor boundary
    else {
        m_factor_min     = min / m_scale;
        m_has_factor_min = true;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set maximum parameter boundary
 *
 * @param[in] max Maximum parameter boundary.
 *
 * @exception GException::invalid_argument
 *            Maximum parameter boundary is smaller than the parameter value
 *
 * Sets the maximum parameter boundary.
 ***************************************************************************/
void GOptimizerPar::max(const double& max)
{
    // Check if maximum is smaller than value
    if (max < value()) {
        std::string msg = "Specified maximum parameter value "+
                          gammalib::str(max)+" is smaller than the parameter "
                          "value "+gammalib::str(value())+".";
        throw GException::invalid_argument(G_MAX, msg);
    }

    // If scale is negative then set minimum factor boundary
    if (m_scale < 0) {
        m_factor_min     = max / m_scale;
        m_has_factor_min = true;
    }

    // ... otherwise set maximum factor boundary
    else {
        m_factor_max     = max / m_scale;
        m_has_factor_max = true;
    }

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
 * Sets the value factor of the parameter. The method makes sure that
 * none of the boundaries is violated. Otherwise, exceptions are thrown.
 ***************************************************************************/
void GOptimizerPar::factor_value(const double& value)
{
    // If there is a minimum boundary and if value is below this boundary
    // then throw an exception
    if (has_factor_min() && (value < factor_min())) {
        std::string msg = "Specified value factor "+gammalib::str(value)+
                          " is smaller than the minimum boundary "+
                          gammalib::str(factor_min())+".";
        throw GException::invalid_argument(G_FACTOR_VALUE, msg);
    }

    // If there is a maximum boundary and if value is above this boundary
    // then throw an exception
    if (has_factor_max() && (value > factor_max())) {
        std::string msg = "Specified value factor "+gammalib::str(value)+
                          " is larger than the maximum boundary "+
                          gammalib::str(factor_max())+".";
        throw GException::invalid_argument(G_FACTOR_VALUE, msg);
    }

    // Assign value
    m_factor_value = value;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set minimum parameter boundary factor
 *
 * @param[in] min Minimum parameter boundary factor.
 *
 * @exception GException::invalid_argument
 *            Minimum parameter boundary factor larger than value factor.
 *
 * Sets the minimum boundary factor of the parameter. The method makes
 * sure that the minimum is not larger than the actual value factor.
 * Otherwise, an exception is thrown.
 ***************************************************************************/
void GOptimizerPar::factor_min(const double& min)
{
    // Check if minimum is larger than value
    if (min > m_factor_value) {
        std::string msg = "Specified minimum factor "+gammalib::str(min)+
                          " is larger than the value factor "+
                          gammalib::str(m_factor_value)+".";
        throw GException::invalid_argument(G_FACTOR_VALUE, msg);
    }

    // Set value and flag
    m_factor_min     = min;
    m_has_factor_min = true;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set maximum parameter boundary factor
 *
 * @param[in] max Maximum parameter boundary factor.
 *
 * @exception GException::invalid_argument
 *            Maximum parameter boundary factor smaller than value factor.
 *
 * Sets the maximum boundary factor of the parameter. The method makes
 * sure that the maximum is not smaller than the actual value factor.
 * Otherwise, an exception is thrown.
 ***************************************************************************/
void GOptimizerPar::factor_max(const double& max)
{
    // Check if maximum is smaller than value
    if (max < m_factor_value) {
        std::string msg = "Specified maximum factor "+gammalib::str(max)+
                          " is smaller than the value factor "+
                          gammalib::str(m_factor_value)+".";
        throw GException::invalid_argument(G_FACTOR_VALUE, msg);
    }

    // Set value and flag
    m_factor_max     = max;
    m_has_factor_max = true;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set minimum and maximum parameter boundary factors
 *
 * @param[in] min Minimum parameter boundary factor.
 * @param[in] max Maximum parameter boundary factor.
 *
 * Sets the minimum and maximum parameter boundary factors. The method calls
 * the factor_min() and factor_max() methods which perform validity checking
 * of the arguments.
 ***************************************************************************/
void GOptimizerPar::factor_range(const double& min, const double& max)
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
 * Sets the scale factor of the parameter. All parameter attributes are
 * rescaled accordingly.
 *
 * An exception is thrown if a scale factor of 0 is specified.
 ***************************************************************************/
void GOptimizerPar::scale(const double& scale)
{
    // Make sure that scale is not zero
    if (scale == 0.0) {
        std::string msg = "Specified parameter scale factor of 0.\n"
                          "Parameters need a non-zero scale factor.";
        throw GException::invalid_argument(G_SCALE, msg);
    }

    // Set rescaling
    double rescale = m_scale/scale;

    // Set new scale factor
    m_scale = scale;

    // Set values, error, gradient, min and max
    m_factor_value    *= rescale;
    m_factor_error    *= rescale;
    m_factor_gradient /= rescale;
    m_factor_min      *= rescale;
    m_factor_max      *= rescale;

    // If re-scaling changes sign then swap minimum and maximum boundary
    if (rescale < 0.0) {

        // Store maximum in swap space
        bool   bswap     = m_has_factor_max;
        double dswap     = m_factor_max;

        // Set maximum
        m_has_factor_max = m_has_factor_min;
        m_factor_max     = m_factor_min;

        // Set minimum
        m_has_factor_min = bswap;
        m_factor_min     = dswap;
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
 * The method does nothing if the actual value factor is zero.
 ***************************************************************************/
void GOptimizerPar::autoscale(void)
{
    // Get value
    double value = this->value();

    // Continue only if the value is non-zero
    if (value != 0.0) {

        // Re-scale to value
        this->scale(value);

    } // endif: value was non-zero

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print parameter information
 *
 * @param[in] chatter Chattiness.
 * @return String with parameter information.
 ***************************************************************************/
std::string GOptimizerPar::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append parameter name
        result.append(gammalib::parformat(" "+name()));

        // Append value
        result.append(gammalib::str(value()));

        // For free parameters, append statistical uncertainty
        if (m_free) {
            result.append(" +/- "+gammalib::str(error()));
        }

        // Append parameter limites if they exist
        if (has_min() && has_max()) {
            result.append(" ["+gammalib::str(min()) + ","+gammalib::str(max())+"]");
        }
        else if (has_min()) {
            result.append(" ["+gammalib::str(min()) + ",infty[");
        }
        else if (has_max()) {
            result.append(" ]-infty,"+gammalib::str(max())+"]");
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
        result.append(",scale="+gammalib::str(m_scale));

        // Signal if parameter has analytic gradient
        if (m_has_grad) {
            result.append(",gradient)");
        }
        else {
            result.append(")");
        }

    } // endif: chatter was not silent

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
void GOptimizerPar::init_members(void)
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
    m_has_factor_min  = false;
    m_has_factor_max  = false;
    m_has_grad        = false;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] par Model parameter.
 ***************************************************************************/
void GOptimizerPar::copy_members(const GOptimizerPar& par)
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
    m_has_factor_min  = par.m_has_factor_min;
    m_has_factor_max  = par.m_has_factor_max;
    m_has_grad        = par.m_has_grad;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GOptimizerPar::free_members(void)
{
    // Return
    return;
}
