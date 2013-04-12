/***************************************************************************
 *       GModelSpatialDiffuseConst.cpp - Spatial isotropic model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GModelSpatialDiffuseConst.cpp
 * @brief Isotropic spatial model class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelSpatialDiffuseConst.hpp"
#include "GModelSpatialRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpatialDiffuseConst g_spatial_const_seed;
const GModelSpatialRegistry     g_spatial_const_registry(&g_spatial_const_seed);

/* __ Method name definitions ____________________________________________ */
#define G_MC                           "GModelSpatialDiffuseConst::mc(GRan&)"
#define G_READ                "GModelSpatialDiffuseConst::read(GXmlElement&)"
#define G_WRITE              "GModelSpatialDiffuseConst::write(GXmlElement&)"

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
GModelSpatialDiffuseConst::GModelSpatialDiffuseConst(void) :
                           GModelSpatialDiffuse()
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
 * Constructs isotropic spatial model by extracting information from an XML
 * element. See the read() method for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelSpatialDiffuseConst::GModelSpatialDiffuseConst(const GXmlElement& xml) :
                           GModelSpatialDiffuse()
{
    // Initialise members
    init_members();

    // Read information from XML element
    read(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Value constructor
 *
 * @param[in] value Isotropic value.
 *
 * Constructs isotropic spatial model by assigning the value of the diffuse
 * emission. This constructor explicitly sets the m_value parameter of the
 * model.
 ***************************************************************************/
GModelSpatialDiffuseConst::GModelSpatialDiffuseConst(const double& value) :
                           GModelSpatialDiffuse()
{
    // Initialise members
    init_members();

    // Set parameter
    m_value.value(value);

    // Perform autoscaling of parameter
    autoscale();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Isotropic spatial model.
 ***************************************************************************/
GModelSpatialDiffuseConst::GModelSpatialDiffuseConst(const GModelSpatialDiffuseConst& model) :
                           GModelSpatialDiffuse(model)
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
GModelSpatialDiffuseConst::~GModelSpatialDiffuseConst(void)
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
 * @return Isotropic spatial model.
 ***************************************************************************/
GModelSpatialDiffuseConst& GModelSpatialDiffuseConst::operator=(const GModelSpatialDiffuseConst& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelSpatialDiffuse::operator=(model);

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
 * @brief Clear isotropic spatial model
 ***************************************************************************/
void GModelSpatialDiffuseConst::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelSpatialDiffuse::free_members();
    this->GModelSpatial::free_members();

    // Initialise members
    this->GModelSpatial::init_members();
    this->GModelSpatialDiffuse::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone isotropic spatial model
 *
 * @return Pointer to deep copy of isotropic spatial model.
 ***************************************************************************/
GModelSpatialDiffuseConst* GModelSpatialDiffuseConst::clone(void) const
{
    // Clone isotropic spatial model
    return new GModelSpatialDiffuseConst(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] photon Incident photon (ignored).
 * @return Model value.
 *
 * Evaluates the spatial part for an isotropic source model. By definition
 * this value is independent from the sky direction.
 ***************************************************************************/
double GModelSpatialDiffuseConst::eval(const GPhoton& photon) const
{
    // Return value
    return (m_value.value());
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] photon Incident photon (ignored).
 * @return Model value.
 *
 * Evaluates the spatial part for an isotropic source model and set the
 * parameter gradient. By definition, the value and gradient is independent
 * from the sky direction. The value is 1, the parameter gradient is 0.
 ***************************************************************************/
double GModelSpatialDiffuseConst::eval_gradients(const GPhoton& photon) const
{
    // Compute function value
    double value = m_value.value();

    // Compute partial derivatives of the parameter values
    double g_norm = (m_value.isfree()) ? m_value.scale() : 0.0;

    // Set gradient (circumvent const correctness)
    const_cast<GModelSpatialDiffuseConst*>(this)->m_value.factor_gradient(g_norm);

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Return MC sky direction
 *
 * @param[in] energy Photon energy (ignored).
 * @param[in] time Photon arrival time (ignored).
 * @param[in,out] ran Random number generator.
 * @return Sky direction.
 *
 * Returns an arbitrary position on the celestial sphere.
 ***************************************************************************/
GSkyDir GModelSpatialDiffuseConst::mc(const GEnergy& energy,
                                      const GTime&   time,
                                      GRan&          ran) const
{
    // Simulate Right Ascension and Declination
    double ra  = gammalib::twopi * ran.uniform();
    double dec = std::acos(1.0 - 2.0 * ran.uniform());

    // Set sky direction
    GSkyDir dir;
    dir.radec(ra, dec);

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
 * element is expected to have the following format:
 *
 *     <spatialModel type="ConstantValue">
 *       <parameter name="Value" scale="1" value="1" min="1"  max="1" free="0"/>
 *     </spatialModel>
 *
 ***************************************************************************/
void GModelSpatialDiffuseConst::read(const GXmlElement& xml)
{
    // Clear model
    clear();

    // Verify that XML element has exactly 1 parameters
    if (xml.elements() != 1 || xml.elements("parameter") != 1) {
        throw GException::model_invalid_parnum(G_READ, xml,
              "Isotropic source model requires exactly 1 parameter.");
    }

    // Get pointer on model parameter
    const GXmlElement* par = xml.element("parameter", 0);

    // Get value
    if (par->attribute("name") == "Value") {
        m_value.read(*par);
    }
    else {
        throw GException::model_invalid_parnames(G_READ, xml,
              "Isotropic source model requires \"Value\" parameter.");
    }

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
 * element will have the following format:
 *
 *     <spatialModel type="ConstantValue">
 *       <parameter name="Value" scale="1" value="1" min="1"  max="1" free="0"/>
 *     </spatialModel>
 *
 ***************************************************************************/
void GModelSpatialDiffuseConst::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "") {
        xml.attribute("type", "ConstantValue");
    }

    // Verify model type
    if (xml.attribute("type") != "ConstantValue") {
        throw GException::model_invalid_spatial(G_WRITE, xml.attribute("type"),
              "Spatial model is not of type \"ConstantValue\".");
    }

    // If XML element has 0 nodes then append 1 parameter node
    if (xml.elements() == 0) {
        xml.append(GXmlElement("parameter name=\"Value\""));
    }

    // Verify that XML element has exactly 1 parameter
    if (xml.elements() != 1 || xml.elements("parameter") != 1) {
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Isotropic source model requires exactly 1 parameter.");
    }

    // Get pointers on both model parameters
    GXmlElement* par = xml.element("parameter", 0);

    // Set or update parameter
    if (par->attribute("name") == "Value") {
        m_value.write(*par);
    }
    else {
        throw GException::model_invalid_parnames(G_WRITE, xml,
              "Isotropic source model requires \"Value\" parameter.");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print isotropic source model information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing model information.
 ***************************************************************************/
std::string GModelSpatialDiffuseConst::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpatialDiffuseConst ===");

        // Append parameters
        result.append("\n"+gammalib::parformat("Number of parameters"));
        result.append(gammalib::str(size()));
        for (int i = 0; i < size(); ++i) {
            result.append("\n"+m_pars[i]->print(chatter));
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
void GModelSpatialDiffuseConst::init_members(void)
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
void GModelSpatialDiffuseConst::copy_members(const GModelSpatialDiffuseConst& model)
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
void GModelSpatialDiffuseConst::free_members(void)
{
    // Return
    return;
}
