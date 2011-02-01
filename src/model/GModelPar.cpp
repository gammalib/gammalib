/***************************************************************************
 *                   GModelPar.cpp  -  Model parameter class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelPar.cpp
 * @brief GModelPar class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GModelPar.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_REAL_VALUE                         "GModelPar::real_value(double&)"
#define G_REAL_ERROR                         "GModelPar::real_error(double&)"
#define G_VALUE                                   "GModelPar::value(double&)"
#define G_MIN                                       "GModelPar::min(double&)"
#define G_MAX                                       "GModelPar::max(double&)"
#define G_RANGE                          "GModelPar::range(double&, double&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                     GModelPar constructors/destructors                  =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GModelPar::GModelPar(void)
{
    // Initialise private members for clean destruction
    init_members();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] par Parameter from which the instance should be built.
 ***************************************************************************/
GModelPar::GModelPar(const GModelPar& par)
{ 
    // Initialise private members for clean destruction
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
 * @param[in] par Parameter which should be assigned.
 ***************************************************************************/
GModelPar& GModelPar::operator= (const GModelPar& par)
{ 
    // Execute only if object is not identical
    if (this != &par) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
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
 * @brief Set real parameter value
 *
 * @param[in] value Real parameter value.
 *
 * @exception GException::model_invalid_parscale
 *            Parameter has a zero scale factor
 *
 * Assigns the real parameter value by dividing by the parameter scale.
 ***************************************************************************/
void GModelPar::real_value(const double& value)
{
    // Throw error if scale is 0
    if (m_scale == 0.0)
        throw GException::model_invalid_parscale(G_REAL_VALUE, m_scale,
              "Zero scale not allowed for model parameter.");
    
    // Compute value
    double val = value / m_scale;

    // Assign value
    this->value(val);
	
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set real parameter error
 *
 * @param[in] error Real parameter error.
 *
 * @exception GException::model_invalid_parscale
 *            Parameter has a zero scale factor
 *
 * Assigns the real parameter error by dividing by the parameter scale.
 ***************************************************************************/
void GModelPar::real_error(const double& error)
{
    // Throw error if scale is 0
    if (m_scale == 0.0)
        throw GException::model_invalid_parscale(G_REAL_ERROR, m_scale,
              "Zero scale not allowed for model parameter.");
    
    // Compute error
    double err = error / m_scale;

    // Assign error
    this->error(err);
	
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set parameter value
 *
 * @param[in] value Parameter value.
 *
 * @exception GException::out_of_range
 *            Parameter value outside valid range
 ***************************************************************************/
void GModelPar::value(const double& value)
{
    // If there is a minimum boundary and if value is below this boundary
    // then throw an error
    if (m_hasmin && value < m_min)
        throw GException::out_of_range(G_VALUE, value, m_min, m_max);

    // If there is a maximum boundary and if value is above this boundary
    // then throw an error
    if (m_hasmax && value > m_max)
        throw GException::out_of_range(G_VALUE, value, m_min, m_max);

    // Assign value
    m_value = value;
	
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set parameter minimum
 *
 * @param[in] min Parameter minimum.
 ***************************************************************************/
void GModelPar::min(const double& min)
{
    // If minimum is above actual value then throw error
    if (m_value < min)
        throw GException::out_of_range(G_MIN, m_value, min, m_max);

    // Assign minimum
    m_min = min;
    
    // Flag that minimum was set
    m_hasmin = true;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set parameter maximum
 *
 * @param[in] max Parameter maximum.
 ***************************************************************************/
void GModelPar::max(const double& max)
{
    // If maximum is below value then throw error
    if (m_value > max)
        throw GException::out_of_range(G_MAX, m_value, m_min, max);

    // Assign maximum
    m_max = max;
    
    // Flag that maximum was set
    m_hasmax = true;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set parameter minimum and maximum
 *
 * @param[in] min Parameter minimum.
 * @param[in] max Parameter maximum.
 ***************************************************************************/
void GModelPar::range(const double& min, const double& max)
{
    // If maximum is below value then throw error
    if (m_value < min || m_value > max || min > max)
        throw GException::out_of_range(G_RANGE, m_value, min, max);

    // Assign range
    m_min = min;
    m_max = max;
    
    // Flag that range was set
    m_hasmin = true;
    m_hasmax = true;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Extract parameter attributes from XML element
 *
 * @param[in] xml XML element containing parameter attributes.
 *
 * The following parameter attributes may exist (they take the default
 * values when they are not found in the XML element):
 * 'value' (default=0.0), 'error' (default=0.0), 'scale' (default=1.0),
 * 'min' (default=INDEF), 'max' (default=INDEF), 'free' (default=0).
 ***************************************************************************/
void GModelPar::read(const GXmlElement& xml)
{
    // Get value
    std::string arg = xml.attribute("value");
    if (arg != "")
        m_value = todouble(arg);
    else
        m_value = 0.0;

    // Get error
    arg = xml.attribute("error");
    if (arg != "")
        m_error = todouble(arg);
    else
        m_error = 0.0;

    // Get scale factor
    arg = xml.attribute("scale");
    if (arg != "")
        m_scale = todouble(arg);
    else
        m_scale = 1.0;

    // Get min
    arg = xml.attribute("min");
    if (arg != "")
        min(todouble(arg));

    // Get max
    arg = xml.attribute("max");
    if (arg != "")
        max(todouble(arg));

    // Get free
    if (xml.attribute("free") == "1" || 
        tolower(xml.attribute("free")) == "true")
        free();
    else
        fix();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set or update parameter attributes in XML element
 *
 * @param[in] xml XML element containing parameter attributes.
 *
 * The following parameter attributes may exist (they take the default
 * values when they are not found in the XML element):
 * 'value' (default=0.0), 'error' (default=0.0), 'scale' (default=1.0),
 * 'min' (default=INDEF), 'max' (default=INDEF), 'free' (default=0).
 ***************************************************************************/
void GModelPar::write(GXmlElement& xml) const
{
    // Set value
    xml.attribute("value", str(m_value));

    // Set error (only if parameter is free)
    if (isfree())
        xml.attribute("error", str(m_error));

    // Set scale
    xml.attribute("scale", str(m_scale));

    // Set minimum
    if (hasmin())
        xml.attribute("min", str(m_min));

    // Set maximum
    if (hasmax())
        xml.attribute("max", str(m_max));

    // Set free/fix flag
    if (isfree())
        xml.attribute("free", "1");
    else
        xml.attribute("free", "0");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print powerlaw information
 ***************************************************************************/
std::string GModelPar::print(void) const
{
    // Initialise result string
    std::string result;

    // Append parameter name
    result.append(parformat(" "+name()));

    // Append real value
    result.append(str(real_value()));

    // For free parameters, append statistical uncertainty
    if (m_free)
        result.append(" +/- "+str(fabs(real_error())));

    // Append parameter limites if they exist
    if (m_hasmin && m_hasmax)
        result.append(" ["+str(real_min()) + ","+str(real_max())+"]");
    else if (m_hasmin)
        result.append(" ["+str(real_min()) + ",infty[");
    else if (m_hasmax)
        result.append(" ]-infty,"+str(real_max())+"]");

    // Append parameter unit
    result.append(" "+m_unit);

    // Signal if parameter was free or fixed
    if (m_free)
        result.append(" (free,");
    else
        result.append(" (fixed,");

    // Append parameter scale
    result.append("scale="+str(m_scale)+")");

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
    m_value    = 0.0;
    m_error    = 0.0;
    m_gradient = 0.0;
    m_min      = 0.0;
    m_max      = 0.0;
    m_scale    = 1.0;
    m_hasmin   = false;
    m_hasmax   = false;
    m_free     = true;
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] par GModelPar members which should be copied.
 ***************************************************************************/
void GModelPar::copy_members(const GModelPar& par)
{
    // Copy members
    m_name    = par.m_name;
    m_unit    = par.m_unit;
    m_value   = par.m_value;
    m_error   = par.m_error;
    m_gradient = par.m_gradient;
    m_min     = par.m_min;
    m_max     = par.m_max;
    m_scale   = par.m_scale;
    m_hasmin  = par.m_hasmin;
    m_hasmax  = par.m_hasmax;
    m_free    = par.m_free;

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


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] par Parameter.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GModelPar& par)
{
     // Write spectrum in output stream
    os << par.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] par Parameter.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GModelPar& par)
{
    // Write spectrum into logger
    log << par.print();

    // Return logger
    return log;
}
