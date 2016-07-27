/***************************************************************************
 *         GModelSpectralConst.cpp - Spectral constant model class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2016 by Juergen Knoedlseder                         *
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
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GRan.hpp"
#include "GModelSpectralConst.hpp"
#include "GModelSpectralRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpectralConst    g_spectral_const_seed1("Constant", "Normalization");
const GModelSpectralRegistry g_spectral_const_registry1(&g_spectral_const_seed1);
#if defined(G_LEGACY_XML_FORMAT)
const GModelSpectralConst    g_spectral_const_seed2("ConstantValue", "Value");
const GModelSpectralRegistry g_spectral_const_registry2(&g_spectral_const_seed2);
#endif

/* __ Method name definitions ____________________________________________ */
#define G_FLUX                "GModelSpectralConst::flux(GEnergy&, GEnergy&)"
#define G_EFLUX              "GModelSpectralConst::eflux(GEnergy&, GEnergy&)"
#define G_MC     "GModelSpectralConst::mc(GEnergy&, GEnergy&, GTime&, GRan&)"
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
 * @brief Model type and parameter name constructor
 *
 * @param[in] type Model type.
 * @param[in] value Name of value parameter.
 ***************************************************************************/
GModelSpectralConst::GModelSpectralConst(const std::string& type,
                                         const std::string& value) :
                     GModelSpectral()
{
    // Initialise members
    init_members();

    // Set model type
    m_type = type;

    // Set parameter names
    m_norm.name(value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs constant spectral model by extracting information from an XML
 * element. See the read() method for more information about the expected
 * structure of the XML element.
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
 * @brief Value constructor
 *
 * @param[in] value Model value (ph/cm2/s/MeV).
 *
 * Constructs constant spectral model by setting the model value.
 ***************************************************************************/
GModelSpectralConst::GModelSpectralConst(const double& value) :
                     GModelSpectral()
{
    // Initialise members
    init_members();

    // Set model value
    m_norm.value(value);

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
 * @param[in] model Constant spectral model.
 * @return Constant spectral model.
 ***************************************************************************/
GModelSpectralConst& GModelSpectralConst::operator=(const GModelSpectralConst& model)
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
 * @brief Clear constant spectral model
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
 * @brief Clone constant spectral model
 *
 * @return Pointer to deep copy of constant spectral model.
 ***************************************************************************/
GModelSpectralConst* GModelSpectralConst::clone(void) const
{
    // Clone constant spectral model
    return new GModelSpectralConst(*this);
}


/***********************************************************************//**
 * @brief Evaluate model value
 *
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @return Model value (ph/cm2/s/MeV).
 *
 * Evaluates
 *
 * \f[
 *    S_{\rm E}(E | t) = {\tt m\_norm}
 * \f]
 *
 * where
 * \f${\tt m\_norm}\f$ is the normalization constant in units of 
 * ph/cm2/s/MeV.
 ***************************************************************************/
double GModelSpectralConst::eval(const GEnergy& srcEng,
                                 const GTime&   srcTime) const
{
    // Compute function value
    double value = m_norm.value();

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate model value and gradient
 *
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @return Model value (ph/cm2/s/MeV).
 *
 * Evaluates
 *
 * \f[
 *    S_{\rm E}(E | t) = {\tt m\_norm}
 * \f]
 *
 * where
 * \f${\tt m\_norm}\f$ is the normalization constant in units of 
 * ph/cm2/s/MeV.
 *
 * The method also evaluates the partial derivative of the model with respect
 * to the m_norm parameter using
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta {\tt m\_norm}} = 1
 * \f]
 ***************************************************************************/
double GModelSpectralConst::eval_gradients(const GEnergy& srcEng,
                                           const GTime&   srcTime)
{
    // Compute function value
    double value = m_norm.value();

    // Compute partial derivatives of the parameter values
    double g_norm = (m_norm.is_free()) ? m_norm.scale() : 0.0;

    // Set factor gradient (the parameter gradient is obtained by dividing
    // the factor gradient by the scale factor)
    m_norm.factor_gradient(g_norm);

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Returns model photon flux between [emin, emax] (ph/cm2/s)
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 * @return Photon flux (ph/cm2/s).
 *
 * Computes
 *
 * \f[
 *    \int_{\tt emin}^{\tt emax} S_{\rm E}(E | t) dE
 * \f]
 *
 * where
 * - [@p emin, @p emax] is an energy interval, and
 * - \f$S_{\rm E}(E | t)\f$ is the spectral model (ph/cm2/s/MeV).
 * The integration is done analytically.
 ***************************************************************************/
double GModelSpectralConst::flux(const GEnergy& emin,
                                 const GEnergy& emax) const
{
    // Initialise flux
    double flux = 0.0;
    
    // Compute only if integration range is valid
    if (emin < emax) {

        // Compute flux for a constant model
        flux = m_norm.value() * (emax.MeV() - emin.MeV());
    
    } // endif: integration range was valid

    // Return
    return flux;
}


/***********************************************************************//**
 * @brief Returns model energy flux between [emin, emax] (erg/cm2/s)
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 * @return Energy flux (erg/cm2/s).
 *
 * Computes
 *
 * \f[
 *    \int_{\tt emin}^{\tt emax} S_{\rm E}(E | t) E \, dE
 * \f]
 *
 * where
 * - [@p emin, @p emax] is an energy interval, and
 * - \f$S_{\rm E}(E | t)\f$ is the spectral model (ph/cm2/s/MeV).
 * The integration is done analytically.
 ***************************************************************************/
double GModelSpectralConst::eflux(const GEnergy& emin,
                                  const GEnergy& emax) const
{
    // Initialise flux
    double eflux = 0.0;
    
    // Compute only if integration range is valid
    if (emin < emax) {

        // Compute flux for a constant model
        eflux = m_norm.value() * 0.5 * (emax.MeV()*emax.MeV() - 
                                        emin.MeV()*emin.MeV());

        // Convert from MeV/cm2/s to erg/cm2/s
        eflux *= gammalib::MeV2erg;
    
    } // endif: integration range was valid

    // Return
    return eflux;
}


/***********************************************************************//**
 * @brief Returns MC energy between [emin, emax]
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 * @param[in] time True photon arrival time.
 * @param[in,out] ran Random number generator.
 * @return Energy.
 *
 * @exception GException::erange_invalid
 *            Energy range is invalid (emin < emax required).
 *
 * Returns Monte Carlo energy by randomly drawing from a constant between
 * the minimum and maximum photon energy.
 ***************************************************************************/
GEnergy GModelSpectralConst::mc(const GEnergy& emin,
                                const GEnergy& emax,
                                const GTime&   time,
                                GRan&          ran) const
{
    // Throw an exception if energy range is invalid
    if (emin >= emax) {
        throw GException::erange_invalid(G_MC, emin.MeV(), emax.MeV(),
              "Minimum energy < maximum energy required.");
    }

    // Allocate energy
    GEnergy energy;

    // Get uniform value between 0 and 1
    double u = ran.uniform();

    // Map into [emin,emax] range
    energy.MeV((emax.MeV() - emin.MeV()) * u + emin.MeV());

    // Return energy
    return energy;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * Reads the spectral information from an XML element.
 ***************************************************************************/
void GModelSpectralConst::read(const GXmlElement& xml)
{
    // Get parameter pointers
    const GXmlElement* norm = gammalib::xml_get_par(G_READ, xml, m_norm.name());

    // Read parameters
    m_norm.read(*norm);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::model_invalid_spectral
 *            Existing XML element is not of type "ConstantValue"
 *
 * Writes the spectral information into an XML element.
 ***************************************************************************/
void GModelSpectralConst::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "") {
        xml.attribute("type", type());
    }

    // Verify model type
    if (xml.attribute("type") != type()) {
        throw GException::model_invalid_spectral(G_WRITE, xml.attribute("type"),
              "Spectral model is not of type \""+type()+"\".");
    }

    // Get normalisation parameter
    GXmlElement* par = gammalib::xml_need_par(G_WRITE, xml, m_norm.name());

    // Write normalisation parameter
    m_norm.write(*par);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print spectral model information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing spectral model information.
 ***************************************************************************/
std::string GModelSpectralConst::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpectralConst ===");

        // Append model content
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
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModelSpectralConst::init_members(void)
{
    // Initialise model type
    m_type = "Constant";

    // Initialise constant normalisation
    m_norm.clear();
    m_norm.name("Normalization");
    m_norm.scale(1.0);
    m_norm.value(1.0);         // default: 1.0
    m_norm.range(0.0, 1000.0); // range:   [0, 1000]
    m_norm.free();
    m_norm.gradient(0.0);
    m_norm.has_grad(true);

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
    m_type = model.m_type;
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
