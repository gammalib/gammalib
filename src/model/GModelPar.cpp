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
GModelPar::GModelPar(void) : GOptimizerPar()
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
GModelPar::GModelPar(const std::string& name, const double& value) :
           GOptimizerPar(name, value)
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
GModelPar::GModelPar(const std::string& name,
                     const double&      factor,
                     const double&      scale) :
           GOptimizerPar(name, factor, scale)
                             
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] par Model parameter.
 ***************************************************************************/
GModelPar::GModelPar(const GModelPar& par) : GOptimizerPar(par)
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

        // Copy base class members
        this->GOptimizerPar::operator=(par);

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(par);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

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
        m_factor_value = gammalib::todouble(arg);
    }
    else {
        m_factor_value = 0.0;
    }

    // Get error
    arg = xml.attribute("error");
    if (arg != "") {
        m_factor_error = gammalib::todouble(arg);
    }
    else {
        m_factor_error = 0.0;
    }

    // Get scale factor
    arg = xml.attribute("scale");
    if (arg != "") {
        m_scale = gammalib::todouble(arg);
    }
    else {
        m_scale = 1.0;
    }

    // Get min
    arg = xml.attribute("min");
    if (arg != "") {
        m_factor_min = gammalib::todouble(arg);
        m_has_min     = true;
    }
    else {
        remove_min();
    }

    // Get max
    arg = xml.attribute("max");
    if (arg != "") {
        m_factor_max = gammalib::todouble(arg);
        m_has_max     = true;
    }
    else {
        remove_max();
    }

    // Get free
    if (xml.attribute("free") == "1" || 
        gammalib::tolower(xml.attribute("free")) == "true") {
        free();
    }
    else {
        fix();
    }

    // If there is a minimum and maximum, make sure that the maximum is
    // not smaller than the minimum
    if (m_has_min && m_has_max) {
        if (m_factor_min > m_factor_max) {
            std::string msg = "The model parameter \""+m_name+
                              "\" in the XML document has a minimum boundary "+
                              gammalib::str(m_factor_min)+
                              " that is larger than the maximum boundary "+
                              gammalib::str(m_factor_max)+".\n"+xml.print();
            throw GException::invalid_value(G_READ, msg);
        }
    }

    // If there is a minimum, make sure that the value is not below it
    if (m_has_min && m_factor_value < m_factor_min) {
        std::string msg = "The model parameter \""+m_name+
                          "\" in the XML document has a value "+
                            gammalib::str(m_factor_value)+
                            " that is smaller than the minimum boundary "+
                            gammalib::str(m_factor_min)+".\n"+xml.print();
        throw GException::invalid_value(G_READ, msg);
    }

    // If there is a maximum, make sure that the value is not above it
    if (m_has_max && m_factor_value > m_factor_max) {
        std::string msg = "The model parameter \""+m_name+
                          "\" in the XML document has a value "+
                            gammalib::str(m_factor_value)+
                            " that is larger than the maximum boundary "+
                            gammalib::str(m_factor_max)+".\n"+xml.print();
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
    xml.attribute("value", gammalib::str(m_factor_value));

    // Set error (only if parameter is free)
    if (is_free()) {
        xml.attribute("error", gammalib::str(m_factor_error));
    }

    // Set scale
    xml.attribute("scale", gammalib::str(m_scale));

    // Set minimum
    if (has_min()) {
        xml.attribute("min", gammalib::str(m_factor_min));
    }

    // Set maximum
    if (has_max()) {
        xml.attribute("max", gammalib::str(m_factor_max));
    }

    // Set free/fix flag
    if (is_free()) {
        xml.attribute("free", "1");
    }
    else {
        xml.attribute("free", "0");
    }

    // Return
    return;
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
