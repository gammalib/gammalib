/***************************************************************************
 *        GModelSpectralConst.cpp  -  Spectral constant model class        *
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
 * @file GModelSpectralConst.cpp
 * @brief Constant spectral model class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelSpectralConst.hpp"
#include "GModelSpectralRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpectralConst    g_spectral_const_seed;
const GModelSpectralRegistry g_spectral_const_registry(&g_spectral_const_seed);

/* __ Method name definitions ____________________________________________ */
#define G_MC             "GModelSpectralConst::mc(GEnergy&, GEnergy&, GRan&)"
#define G_READ                      "GModelSpectralConst::read(GXmlElement&)"
#define G_WRITE                    "GModelSpectralConst::write(GXmlElement&)"

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
GModelSpectralConst::GModelSpectralConst(void) : GModelSpectral()
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
 * Creates instance of a constant spectral model by extracting information
 * from an XML element. See GModelSpectralConst::read() for more information
 * about the expected structure of the XML element.
 ***************************************************************************/
GModelSpectralConst::GModelSpectralConst(const GXmlElement& xml) :
                     GModelSpectral()
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
 * @param[in] model Spectral constant model.
 ***************************************************************************/
GModelSpectralConst::GModelSpectralConst(const GModelSpectralConst& model) :
                     GModelSpectral(model)
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
GModelSpectralConst::~GModelSpectralConst(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] model Spectral constant model.
 ***************************************************************************/
GModelSpectralConst& GModelSpectralConst::operator= (const GModelSpectralConst& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelSpectral::operator=(model);

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
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
***************************************************************************/
void GModelSpectralConst::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelSpectral::free_members();

    // Initialise members
    this->GModelSpectral::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GModelSpectralConst* GModelSpectralConst::clone(void) const
{
    return new GModelSpectralConst(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] srcEng True energy of photon.
 *
 * The spectral model is defined as
 * \f[I(E)=norm\f]
 * where
 * \f$norm\f$ is the normalization of the function.
 ***************************************************************************/
double GModelSpectralConst::eval(const GEnergy& srcEng) const
{
    // Compute function value
    double value = norm();

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] srcEng True energy of photon.
 *
 * The spectral model is defined as
 * \f[I(E)=norm\f]
 * where
 * \f$norm=n_s n_v\f$ is the normalization of the function.
 * Note that the normalization is factorised into a scaling factor and a
 * value and that the method is expected to return the gradient with respect
 * to the parameter value \f$n_v\f$.
 *
 * The partial derivative of the normalization value is given by
 * \f[dI/dn_v=n_s\f]
 ***************************************************************************/
double GModelSpectralConst::eval_gradients(const GEnergy& srcEng) const
{
    // Compute function value
    double value = norm();

    // Compute partial derivatives of the parameter values
    double g_norm = (m_norm.isfree()) ? m_norm.scale() : 0.0;

    // Set gradients (circumvent const correctness)
    ((GModelSpectralConst*)this)->m_norm.gradient(g_norm);

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Returns model photon flux between [emin, emax] (units: ph/cm2/s)
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 *
 * Computes
 * \f[\int_{E_{\rm min}}^{E_{\rm max}} I(E) dE\f]
 * where
 * \f$E_{\rm min}\f$ and \f$E_{\rm max}\f$ are the minimum and maximum
 * energy, respectively, and
 * \f$I(E)\f$ is the spectral model (units: ph/cm2/s/MeV).
 * The integration is done analytically.
 ***************************************************************************/
double GModelSpectralConst::flux(const GEnergy& emin, const GEnergy& emax) const
{
    // Compute flux for a constant model
    double flux = norm() * (emax.MeV() - emin.MeV());

    // Return
    return flux;
}


/***********************************************************************//**
 * @brief Returns model energy flux between [emin, emax] (units: erg/cm2/s)
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 *
 * Computes
 * \f[\int_{E_{\rm min}}^{E_{\rm max}} I(E) E dE\f]
 * where
 * \f$E_{\rm min}\f$ and \f$E_{\rm max}\f$ are the minimum and maximum
 * energy, respectively, and
 * \f$I(E)\f$ is the spectral model (units: ph/cm2/s/MeV).
 * The integration is done analytically.
 ***************************************************************************/
double GModelSpectralConst::eflux(const GEnergy& emin, const GEnergy& emax) const
{
    // Compute flux for a constant model
    double flux = norm() * 0.5 * (emax.MeV()*emax.MeV() - emin.MeV()*emin.MeV());

    // Convert from MeV/cm2/s to erg/cm2/s
    flux *= MeV2erg;

    // Return
    return flux;
}


/***********************************************************************//**
 * @brief Returns MC energy between [emin, emax]
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 * @param[in] ran Random number generator.
 *
 * @todo To be implemented
 ***************************************************************************/
GEnergy GModelSpectralConst::mc(const GEnergy& emin, const GEnergy& emax,
                                GRan& ran) const
{
    // Allocate energy
    GEnergy energy;

    // Dump warning that method is not yet implemented
    throw GException::feature_not_implemented(G_MC);

    // Return energy
    return energy;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element containing power law model information.
 *
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter name found in XML element.
 *
 * Read the function information from an XML element and load the nodes
 * from the associated file. The XML element is required to have a parameter
 * named "Normalization" or "Value".
 ***************************************************************************/
void GModelSpectralConst::read(const GXmlElement& xml)
{
    // Verify that XML element has exactly 1 parameter
    if (xml.elements() != 1 || xml.elements("parameter") != 1)
        throw GException::model_invalid_parnum(G_READ, xml,
              "Spectral constant requires exactly 1 parameter.");

    // Get parameter element
    GXmlElement* par = (GXmlElement*)xml.element("parameter", 0);

    // Get value
    if (par->attribute("name") == "Normalization" ||
        par->attribute("name") == "Value")
        m_norm.read(*par);
    else
        throw GException::model_invalid_parnames(G_READ, xml,
                          "Spectral constant requires either"
                          " \"Normalization\" or \"Value\" parameter.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * @exception GException::model_invalid_spectral
 *            Existing XML element is not of type "ConstantValue"
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters or nodes found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Write the spectral constant information into an XML element. The XML
 * element has to be of type "ConstantValue" and will have 1 parameter leaf
 * named "Normalization" or "Value" (default).
 ***************************************************************************/
void GModelSpectralConst::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "")
        xml.attribute("type", "ConstantValue");

    // Verify model type
    if (xml.attribute("type") != "ConstantValue")
        throw GException::model_invalid_spectral(G_WRITE, xml.attribute("type"),
              "Spectral model is not of type \"ConstantValue\".");

    // If XML element has 0 nodes then append parameter node. The name
    // of the node is "Value" as this is the Fermi-LAT standard.
    // We thus assure that the XML files will be compatible with
    // Fermi-LAT.
    if (xml.elements() == 0)
        xml.append(new GXmlElement("parameter name=\"Value\""));

    // Verify that XML element has exactly 1 parameter
    if (xml.elements() != 1 || xml.elements("parameter") != 1)
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Spectral constant requires exactly 1 parameter.");

    // Get parameter element
    GXmlElement* par = (GXmlElement*)xml.element("parameter", 0);

    // Set parameyter
    if (par->attribute("name") == "Normalization" ||
        par->attribute("name") == "Value")
        m_norm.write(*par);
    else
        throw GException::model_invalid_parnames(G_WRITE, xml,
                          "Spectral constant requires either"
                          " \"Normalization\" or \"Value\" parameter.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print powerlaw information
 ***************************************************************************/
std::string GModelSpectralConst::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GModelSpectralConst ===\n");
    result.append(parformat("Number of parameters")+str(size()));
    for (int i = 0; i < size(); ++i)
        result.append("\n"+m_pars[i]->print());

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModelSpectralConst::init_members(void)
{
    // Initialise constant normalisation
    m_norm.clear();
    m_norm.name("Value");
    m_norm.scale(1.0);
    m_norm.value(1.0);         // default: 1.0
    m_norm.range(0.0, 1000.0); // range:   [0, 1000]
    m_norm.free();
    m_norm.gradient(0.0);
    m_norm.hasgrad(true);

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Spectral constant model.
 ***************************************************************************/
void GModelSpectralConst::copy_members(const GModelSpectralConst& model)
{
    // Copy members
    m_norm = model.m_norm;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralConst::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
