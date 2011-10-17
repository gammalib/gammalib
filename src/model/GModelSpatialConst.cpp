/***************************************************************************
 *         GModelSpatialConst.cpp  -  Spatial isotropic model class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @file GModelSpatialConst.cpp
 * @brief Isotropic spatial model class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelSpatialConst.hpp"
#include "GModelSpatialRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpatialConst    g_spatial_const_seed;
const GModelSpatialRegistry g_spatial_const_registry(&g_spatial_const_seed);

/* __ Method name definitions ____________________________________________ */
#define G_MC                                  "GModelSpatialConst::mc(GRan&)"
#define G_READ                       "GModelSpatialConst::read(GXmlElement&)"
#define G_WRITE                     "GModelSpatialConst::write(GXmlElement&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GModelSpatialConst::GModelSpatialConst(void) : GModelSpatial()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Creates instance of isotropic spatial model by extracting information from
 * an XML element. See GModelSpatialConst::read() for more information about
 * the expected structure of the XML element.
 ***************************************************************************/
GModelSpatialConst::GModelSpatialConst(const GXmlElement& xml) :
                    GModelSpatial()
{
    // Initialise members
    init_members();

    // Read information from XML element
    read(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Isotropic spatial model.
 ***************************************************************************/
GModelSpatialConst::GModelSpatialConst(const GModelSpatialConst& model) :
                    GModelSpatial(model)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(model);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GModelSpatialConst::~GModelSpatialConst(void)
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
 * @param[in] model Isotropic spatial model.
 ***************************************************************************/
GModelSpatialConst& GModelSpatialConst::operator= (const GModelSpatialConst& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelSpatial::operator=(model);

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(model);

    } // endif: object was not identical

    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                            Public methods                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
***************************************************************************/
void GModelSpatialConst::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelSpatial::free_members();

    // Initialise members
    this->GModelSpatial::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GModelSpatialConst* GModelSpatialConst::clone(void) const
{
    return new GModelSpatialConst(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] srcDir True photon arrival direction.
 *
 * Evaluates the spatial part for an isotropic source model. By definition
 * this value is independent from the sky direction and is unity.
 ***************************************************************************/
double GModelSpatialConst::eval(const GSkyDir& srcDir) const
{
    // Return value
    return 1.0;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] srcDir True photon arrival direction.
 *
 * Evaluates the spatial part for an isotropic source model and set the
 * parameter gradient. By definition, the value and gradient is independent
 * from the sky direction. The value is 1, the parameter gradient is 0.
 ***************************************************************************/
double GModelSpatialConst::eval_gradients(const GSkyDir& srcDir) const
{
    // Set gradient to 0 (circumvent const correctness)
    ((GModelSpatialConst*)this)->m_value.gradient(0.0);

    // Return value
    return 1.0;
}


/***********************************************************************//**
 * @brief Returns MC sky direction
 *
 * @param[in] ran Random number generator.
 *
 * @exception GException::feature_not_implemented
 *            Method not yet implemented
 *
 * @todo Implement method
 ***************************************************************************/
GSkyDir GModelSpatialConst::mc(GRan& ran) const
{
    // Allocate sky direction
    GSkyDir dir;

    // Dump warning that method is not yet implemented
    throw GException::feature_not_implemented(G_MC);
    
    // Return sky direction
    return dir;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Read the isotropic source model information from an XML element. The XML
 * element is required to have 1 parameter named "Value".
 ***************************************************************************/
void GModelSpatialConst::read(const GXmlElement& xml)
{
    // Verify that XML element has exactly 1 parameters
    if (xml.elements() != 1 || xml.elements("parameter") != 1)
        throw GException::model_invalid_parnum(G_READ, xml,
              "Isotropic source model requires exactly 1 parameter.");

    // Get pointer on model parameter
    GXmlElement* par = (GXmlElement*)xml.element("parameter", 0);

    // Get value
    if (par->attribute("name") == "Value")
        m_value.read(*par);
    else
        throw GException::model_invalid_parnames(G_READ, xml,
              "Isotropic source model requires \"Value\" parameter.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::model_invalid_spatial
 *            Existing XML element is not of type "ConstantValue"
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Write the isotropic source model information into an XML element. The XML
 * element has to be of type "ConstantValue" and will have 1 parameter leaf
 * named "Value".
 ***************************************************************************/
void GModelSpatialConst::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "")
        xml.attribute("type", "ConstantValue");

    // Verify model type
    if (xml.attribute("type") != "ConstantValue")
        throw GException::model_invalid_spatial(G_WRITE, xml.attribute("type"),
              "Spatial model is not of type \"ConstantValue\".");

    // If XML element has 0 nodes then append 1 parameter node
    if (xml.elements() == 0) {
        xml.append(new GXmlElement("parameter name=\"Value\""));
    }

    // Verify that XML element has exactly 1 parameter
    if (xml.elements() != 1 || xml.elements("parameter") != 1)
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Isotropic source model requires exactly 1 parameter.");

    // Get pointers on both model parameters
    GXmlElement* par = (GXmlElement*)xml.element("parameter", 0);

    // Set or update parameter
    if (par->attribute("name") == "Value")
        m_value.write(*par);
    else
        throw GException::model_invalid_parnames(G_WRITE, xml,
              "Isotropic source model requires \"Value\" parameter.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print isotropic source model information
 ***************************************************************************/
std::string GModelSpatialConst::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GModelSpatialConst ===\n");
    result.append(parformat("Number of parameters")+str(size()));
    for (int i = 0; i < size(); ++i)
        result.append("\n"+m_pars[i]->print());

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
void GModelSpatialConst::init_members(void)
{
    // Initialise Value
    m_value.clear();
    m_value.name("Value");
    m_value.fix();
    m_value.value(1.0);
    m_value.scale(1.0);
    m_value.range(0.0, 10.0);
    m_value.gradient(0.0);
    m_value.hasgrad(true);

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Isotropic spatial model.
 ***************************************************************************/
void GModelSpatialConst::copy_members(const GModelSpatialConst& model)
{
    // Copy members
    m_value = model.m_value;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatialConst::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/
