/***************************************************************************
 *         GModelSpectralBrokenPlaw.cpp - Spectral Broken power law model class         *
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
 * @file GModelSpectralBrokenPlaw.cpp
 * @brief Power law spectral broken model class implementation
 * @author Anneli Schulz
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelSpectralBrokenPlaw.hpp"
#include "GModelSpectralRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpectralBrokenPlaw     g_spectral_plaw_seed;
const GModelSpectralRegistry g_spectral_plaw_registry(&g_spectral_plaw_seed);

/* __ Method name definitions ____________________________________________ */
#define G_MC      "GModelSpectralBrokenPlaw::mc(GEnergy&, GEnergy&, GTime&, GRan&)"
#define G_READ                       "GModelSpectralBrokenPlaw::read(GXmlElement&)"
#define G_WRITE                     "GModelSpectralBrokenPlaw::write(GXmlElement&)"

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
GModelSpectralBrokenPlaw::GModelSpectralBrokenPlaw(void) : GModelSpectral()
{
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Parameter constructor
 *
 * @param[in] prefactor Power law pre factor (ph/cm2/s/MeV).
 * @param[in] index1 Power law index1.
 * @param[in] breakenergy break energy.
 * @param[in] index2 Power law index1.
 *
 * Constructs a spectral power law using the model parameters
 * - power law @p prefactor (ph/cm2/s/MeV)
 * - spectral @p index1
 * - breakenergy @p breakenergy
 * - spectral @p index2
 ***************************************************************************/
GModelSpectralBrokenPlaw::GModelSpectralBrokenPlaw(const double&  prefactor,
                                       const double&  index1,
                                       const GEnergy& breakenergy,
                                       const double&  index2) :
                    GModelSpectral()
{
    // Initialise members
    init_members();

    // Set parameters
    m_norm.value(prefactor);
    m_index1.value(index1);
    m_breakenergy.value(breakenergy.MeV());  // Internally stored in MeV
    m_index2.value(index2);

    // Perform autoscaling of parameter
    autoscale();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element containing position information.
 *
 * Constructs a power law spectral model by extracting information from an
 * XML element. See the read() method for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelSpectralBrokenPlaw::GModelSpectralBrokenPlaw(const GXmlElement& xml) : GModelSpectral()
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
 * @param[in] model spectral broken power law model.
 ***************************************************************************/
GModelSpectralBrokenPlaw::GModelSpectralBrokenPlaw(const GModelSpectralBrokenPlaw& model)
                                       : GModelSpectral(model)
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
GModelSpectralBrokenPlaw::~GModelSpectralBrokenPlaw(void)
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
 * @param[in] model spectral broken power law model.
 * @return spectral broken power law model.
 ***************************************************************************/
GModelSpectralBrokenPlaw& GModelSpectralBrokenPlaw::operator=(const GModelSpectralBrokenPlaw& model)
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
 * @brief Clear spectral broken power law model
 ***************************************************************************/
void GModelSpectralBrokenPlaw::clear(void)
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
 * @brief Clone spectral broken power law model
 *
 * @return Pointer to deep copy of spectral broken power law model.
 ***************************************************************************/
GModelSpectralBrokenPlaw* GModelSpectralBrokenPlaw::clone(void) const
{
    // Clone spectral broken power law model
    return new GModelSpectralBrokenPlaw(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @return Model value (ph/cm2/s/MeV).
 *
 * Evaluates
 *
 * \f[
 *    S_{\rm E}(E | t) = {\tt m\_norm}
 *    \left( \frac{E}{\tt m\_breakenergy} \right)^{\tt m\_index}
 * \f]
 *
 * where
 * - \f${\tt m\_norm}\f$ is the normalization or prefactor,
 * - \f${\tt m\_index}\f$ is the spectral index, and
 * - \f${\tt m\_breakenergy}\f$ is the breakenergy energy.
 *
 * @todo The method expects that energy!=0. Otherwise Inf or NaN may result.
 ***************************************************************************/
double GModelSpectralBrokenPlaw::eval(const GEnergy& srcEng,
                                const GTime&   srcTime) const
{
    // Update the evaluation cache
    update_eval_cache(srcEng);

    // Compute function value
    double value = m_norm.value() * m_last_power;

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::isnotanumber(value) || gammalib::isinfinite(value)) {
        std::cout << "*** ERROR: GModelSpectralBrokenPlaw::eval";
        std::cout << "(srcEng=" << srcEng;
        std::cout << ", srcTime=" << srcTime << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", m_norm=" << m_norm.value();
        std::cout << ", m_index1=" << m_index1.value();
        std::cout << ", m_breakenergy=" << m_breakenergy.value();
        std::cout << ", m_index2=" << m_index2.value();
        std::cout << ", m_last_power=" << m_last_power;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @return Model value (ph/cm2/s/MeV).
 *
 * Evaluates
 *
 * \f[
 *    S_{\rm E}(E | t) = {\tt m\_norm}
 *    \left( \frac{E}{\tt m\_breakenergy} \right)^{\tt m\_index}
 * \f]
 *
 * where
 * - \f${\tt m\_norm}\f$ is the normalization or prefactor,
 * - \f${\tt m\_index}\f$ is the spectral index, and
 * - \f${\tt m\_breakenergy}\f$ is the breakenergy energy.
 *
 * The method also evaluates the partial derivatives of the model with
 * respect to the parameters using
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta {\tt m\_norm}} =
 *      \frac{S_{\rm E}(E | t)}{{\tt m\_norm}}
 * \f]
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta {\tt m\_index}} =
 *      S_{\rm E}(E | t) \, \ln(E/{\tt m_breakenergy})
 * \f]
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta {\tt m\_breakenergy}} =
 *      -S_{\rm E}(E | t) \,
 *      \left( \frac{{\tt m\_index}}{{\tt m\_breakenergy}} \right)
 * \f]
 *
 * @todo The method expects that energy!=0. Otherwise Inf or NaN may result.
 ***************************************************************************/
double GModelSpectralBrokenPlaw::eval_gradients(const GEnergy& srcEng,
                                          const GTime&   srcTime)
{
    // Update the evaluation cache
    update_eval_cache(srcEng);

    // Compute function value
    double value = m_norm.value() * m_last_power;

    // Compute partial derivatives of the parameter values
    double g_norm  = (m_norm.isfree())
                     ? m_norm.scale() * m_last_power : 0.0;
    double g_index1 = (m_index1.isfree())
                     ? value * m_index1.scale() * m_last_log_e_norm : 0.0;
    double g_index2 = (m_index2.isfree())
                     ? value * m_index2.scale() * m_last_log_e_norm : 0.0;
    double index_current = (srcEng > m_breakenergy)
                     ? m_last_index1 : m_last_index2;
    double g_breakenergy = (m_breakenergy.isfree())
                     ? -value * index_current / m_breakenergy.factor_value() : 0.0;


    // Set gradients
    m_norm.factor_gradient(g_norm);
    m_index1.factor_gradient(g_index1);
    m_breakenergy.factor_gradient(g_breakenergy);
    m_index2.factor_gradient(g_index2);

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::isnotanumber(value) || gammalib::isinfinite(value)) {
        std::cout << "*** ERROR: GModelSpectralBrokenPlaw::eval_gradients";
        std::cout << "(srcEng=" << srcEng;
        std::cout << ", srcTime=" << srcTime << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", e_norm=" << m_last_e_norm;
        std::cout << ", power=" << m_last_power;
        std::cout << ", g_norm=" << g_norm;
        std::cout << ", g_index1=" << g_index1;
        std::cout << ", g_breakenergy=" << g_breakenergy;
        std::cout << ", g_index2=" << g_index2;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Returns model photon flux between [emin, emax] (units: ph/cm2/s)
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
double GModelSpectralBrokenPlaw::flux(const GEnergy& emin,
                                const GEnergy& emax) const
{
    // Initialise flux
    double flux = 0.0;

    // Compute only if integration range is valid
    if (emin < emax) {
        // first case: emin < breakenergy < emax
        if (emin < m_breakenergy.value() && m_breakenergy.value() < emax) {
            // Compute photon flux
            flux = m_norm.value() *
                   (gammalib::plaw_photon_flux(emin.MeV(),
                                              m_breakenergy.value(),
                                              m_breakenergy.value(),
                                              m_index1.value()) +
                   gammalib::plaw_photon_flux(m_breakenergy.value(),
                                              emax.MeV(),
                                              m_breakenergy.value(),
                                              m_index2.value())
                                              );
        }
        //second case: breakenergy > emax:
        else if (m_breakenergy.value() > emax) {
            flux = m_norm.value() *
                   gammalib::plaw_photon_flux(emin.MeV(),
                                              emax.MeV(),
                                              m_breakenergy.value(),
                                              m_index1.value());
        }
        //third case breakenergy < emin:
        else if (m_breakenergy.value() < emin) {
            flux = m_norm.value() *
                   gammalib::plaw_photon_flux(emin.MeV(),
                                              emax.MeV(),
                                              m_breakenergy.value(),
                                              m_index2.value());
        }
    } // endif: integration range was valid
    // Return flux
    return flux;
}


/***********************************************************************//**
 * @brief Returns model energy flux between [emin, emax] (units: erg/cm2/s)
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
double GModelSpectralBrokenPlaw::eflux(const GEnergy& emin,
                                 const GEnergy& emax) const
{
    // Initialise flux
    double eflux = 0.0;

    // Compute only if integration range is valid
    if (emin < emax) {
        // first case: emin < breakenergy < emax
        if (emin < m_breakenergy.value() && m_breakenergy.value() < emax) {
            eflux = m_norm.value() *
                   (gammalib::plaw_energy_flux(emin.MeV(),
                                              m_breakenergy.value(),
                                              m_breakenergy.value(),
                                              m_index1.value()) +
                   gammalib::plaw_energy_flux(m_breakenergy.value(),
                                              emax.MeV(),
                                              m_breakenergy.value(),
                                              m_index2.value())
                                              );
        }
        //second case: breakenergy > emax:
        else if (m_breakenergy.value() > emax) {
            eflux = m_norm.value() *
                   gammalib::plaw_energy_flux(emin.MeV(),
                                              emax.MeV(),
                                              m_breakenergy.value(),
                                              m_index1.value());
        }
        //third case breakenergy < emin:
        else if (m_breakenergy.value() < emin) {
            eflux = m_norm.value() *
                   gammalib::plaw_energy_flux(emin.MeV(),
                                              emax.MeV(),
                                              m_breakenergy.value(),
                                              m_index2.value());
        }
        eflux *= gammalib::MeV2erg;

    } // endif: integration range was valid

    // Return flux
    return eflux;
}


/***********************************************************************//**
 * @brief Returns Monte Carlo energy between [emin, emax]
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
 * Returns Monte Carlo energy by randomly drawing from a broken power law.
 ***************************************************************************/
GEnergy GModelSpectralBrokenPlaw::mc(const GEnergy& emin,
                               const GEnergy& emax,
                               const GTime&   time,
                               GRan&          ran) const
{
    // Throw an exception if energy range is invalid
    if (emin >= emax) {
        throw GException::erange_invalid(G_MC, emin.MeV(), emax.MeV(),
              "Minimum energy < maximum energy required.");
    }

    // Update cache
    update_mc_cache(emin, emax);

    // Get uniform random number
    double u = ran.uniform();

    // Initialise energy
    double eng;

    // Case A: Index is not -1
    if (index1() != -1.0) {
        if (u > 0.0) {
            eng = std::exp(std::log(u * m_mc_pow_ewidth + m_mc_pow_emin) /
                           m_mc_exponent1);
        }
        else {
            eng = 0.0;
        }
    }

    // Case B: Index is -1
    else {
        eng = std::exp(u * m_mc_pow_ewidth + m_mc_pow_emin);
    }

    // Set energy
    GEnergy energy;
    energy.MeV(eng);

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
 *            Invalid model parameter names found in XML element.
 *
 * Reads the spectral information from an XML element. The format of the XML
 * elements is
 *
 *    <spectrum type="BrokenPowerLaw">
 *      <parameter name="Prefactor" scale=".." value=".." max=".." min=".." free=".."/>
 *      <parameter name="Index1" scale=".." value=".." max=".." min=".." free=".."/>
 *      <parameter name="BreakValue" scale=".." value=".." max=".." min=".." free=".."/>
 *      <parameter name="Index2" scale=".." value=".." max=".." min=".." free=".."/>
 *    </spectrum>
 * @todo Add parameter validity check
 ***************************************************************************/
void GModelSpectralBrokenPlaw::read(const GXmlElement& xml)
{
    // Verify that XML element has exactly 4 parameters
    if (xml.elements() != 4 || xml.elements("parameter") != 4) {
        throw GException::model_invalid_parnum(G_READ, xml,
              "Broken Power law model requires exactly 4 parameters.");
    }

    // Extract model parameters
    int npar[] = {0, 0, 0, 0};
    for (int i = 0; i < 4; ++i) {

        // Get parameter element
        const GXmlElement* par = xml.element("parameter", i);

        // Handle prefactor
        if (par->attribute("name") == "Prefactor") {
            m_norm.read(*par);
            npar[0]++;
        }

        // Handle index1
        else if (par->attribute("name") == "Index1") {
            m_index1.read(*par);
            npar[1]++;
        }

        // Handle breakenergy energy
        else if (par->attribute("name") == "BreakEnergy") {
            m_breakenergy.read(*par);
            npar[2]++;
        }
        // Handle index2
        else if (par->attribute("name") == "Index2") {
            m_index2.read(*par);
            npar[3]++;
                }


    } // endfor: looped over all parameters

    // Verify that all parameters were found
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1 || npar[3] !=1) {
        throw GException::model_invalid_parnames(G_READ, xml,
              "Require \"Prefactor\", \"Index1\", \"Break Energy\" and  \"Index2\" parameters.");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * @exception GException::model_invalid_
 *            Existing XML element is not of type "BrokenPowerLaw"
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters or nodes found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Writes the spectral information into an XML element. The format of the XML
 * element is
 * *    <spectrum type="BrokenPowerLaw">
 *      <parameter name="Prefactor" scale=".." value=".." max=".." min=".." free=".."/>
 *      <parameter name="Index1" scale=".." value=".." max=".." min=".." free=".."/>
 *      <parameter name="BreakValue" scale=".." value=".." max=".." min=".." free=".."/>
 *      <parameter name="Index2" scale=".." value=".." max=".." min=".." free=".."/>
 *    </spectrum>
 ***************************************************************************/
void GModelSpectralBrokenPlaw::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "") {
        xml.attribute("type", "BrokenPowerLaw");
    }

    // Verify model type
    if (xml.attribute("type") != "BrokenPowerLaw") {
        throw GException::model_invalid_spectral(G_WRITE, xml.attribute("type"),
              "Spectral model is not of type \"BrokenPowerLaw\".");
    }

    // If XML element has 0 nodes then append 4 parameter nodes
    if (xml.elements() == 0) {
        xml.append(GXmlElement("parameter name=\"Prefactor\""));
        xml.append(GXmlElement("parameter name=\"Index1\""));
        xml.append(GXmlElement("parameter name=\"BreakEnergy\""));
        xml.append(GXmlElement("parameter name=\"Index2\""));
    }

    // Verify that XML element has exactly 4 parameters
    if (xml.elements() != 4 || xml.elements("parameter") != 4) {
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Broken Power law model requires exactly 4 parameters.");
    }

    // Set or update model parameter attributes
    int npar[] = {0, 0, 0, 0};
    for (int i = 0; i < 4; ++i) {
        // Get parameter element
        GXmlElement* par = xml.element("parameter", i);

        // Handle prefactor
        if (par->attribute("name") == "Prefactor") {
            npar[0]++;
            m_norm.write(*par);
        }
        // Handle index1
        else if (par->attribute("name") == "Index1") {
            npar[1]++;
            m_index1.write(*par);
        }
        // Handle breakenergy energy
        else if (par->attribute("name") == "BreakEnergy") {
            npar[2]++;
            m_breakenergy.write(*par);
        }
        // Handle index2
        else if (par->attribute("name") == "Index2") {
            npar[3]++;
            m_index2.write(*par);
        }
    } // endfor: looped over all parameters

    // Check of all required parameters are present
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1 || npar[3] != 1) {
        throw GException::model_invalid_parnames(G_WRITE, xml,
              "Require \"Prefactor\", \"Index1\", \"Scale\"  and \"Index2\"  parameters.");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print brokenpowerlaw information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing model information.
 ***************************************************************************/
std::string GModelSpectralBrokenPlaw::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpectralBrokenPlaw ===");

        // Append information
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
void GModelSpectralBrokenPlaw::init_members(void)
{
    // Initialise powerlaw normalisation
    m_norm.clear();
    m_norm.name("Prefactor");
    m_norm.unit("ph/cm2/s/MeV");
    m_norm.scale(1.0);
    m_norm.value(1.0);          // default: 1.0
    m_norm.min(0.0);            // min:     0.0
    m_norm.free();
    m_norm.gradient(0.0);
    m_norm.hasgrad(true);

    // Initialise powerlaw index1
    m_index1.clear();
    m_index1.name("Index1");
    m_index1.scale(1.0);
    m_index1.value(-2.0);        // default: -2.0
    m_index1.range(-10.0,+10.0); // range:   [-10,+10]
    m_index1.free();
    m_index1.gradient(0.0);
    m_index1.hasgrad(true);

    // Initialise breakenergy energy
    m_breakenergy.clear();
    m_breakenergy.name("BreakEnergy");
    m_breakenergy.unit("MeV");
    m_breakenergy.scale(1.0);
    m_breakenergy.value(100.0);       // default: 100
    m_breakenergy.fix();
    m_breakenergy.gradient(0.0);
    m_breakenergy.hasgrad(true);

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);
    m_pars.push_back(&m_index1);
    m_pars.push_back(&m_breakenergy);
    m_pars.push_back(&m_index2);

    // Initialise eval cache
    m_last_energy.clear();
    m_last_index1      = 1.0e30;
    m_last_index2      = 1.0e30;
    m_last_breakenergy      = 1.0e30;
    m_last_e_norm     = 0.0;
    m_last_log_e_norm = 0.0;
    m_last_power      = 0.0;

    // Initialise MC cache
    m_mc_emin       = 0.0;
    m_mc_emax       = 0.0;
    m_mc_exponent1   = 0.0;
    m_mc_exponent2   = 0.0;
    m_mc_pow_emin   = 0.0;
    m_mc_pow_ewidth = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model GModelSpectralBrokenPlaw members which should be copied.
 ***************************************************************************/
void GModelSpectralBrokenPlaw::copy_members(const GModelSpectralBrokenPlaw& model)
{
    // Copy members
    m_norm  = model.m_norm;
    m_index1 = model.m_index1;
    m_breakenergy = model.m_breakenergy;
    m_index2 = model.m_index2;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);
    m_pars.push_back(&m_index1);
    m_pars.push_back(&m_breakenergy);
    m_pars.push_back(&m_index2);

    // Copy eval cache
    m_last_energy     = model.m_last_energy;
    m_last_index1      = model.m_last_index1;
    m_last_breakenergy      = model.m_last_breakenergy;
    m_last_index2      = model.m_last_index2;
    m_last_e_norm     = model.m_last_e_norm;
    m_last_log_e_norm = model.m_last_log_e_norm;
    m_last_power      = model.m_last_power;

    // Copy MC cache
    m_mc_emin       = model.m_mc_emin;
    m_mc_emax       = model.m_mc_emax;
    m_mc_exponent1   = model.m_mc_exponent1;
    m_mc_exponent2   = model.m_mc_exponent2;
    m_mc_pow_emin   = model.m_mc_pow_emin;
    m_mc_pow_ewidth = model.m_mc_pow_ewidth;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralBrokenPlaw::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update eval precomputation cache
 *
 * @param[in] energy Energy.
 *
 * Updates the precomputation cache for eval() and eval_gradients() methods.
 ***************************************************************************/
void GModelSpectralBrokenPlaw::update_eval_cache(const GEnergy& energy) const
{
    // Get parameter values (takes 2 multiplications which are difficult
    // to avoid)
    double index1 = m_index1.value();
    double breakenergy = m_breakenergy.value();
    double index2 = m_index2.value();

    // If the energy or one of the parameters index1, index2 or breakenergy energy has
    // changed then recompute the cache
    if ((m_last_energy != energy) ||
        (m_last_index1  != index1)  ||
        (m_last_index2  != index2)  ||
        (m_last_breakenergy  != breakenergy)) {

        // Store actual energy and parameter values
        m_last_energy = energy;
        m_last_index1  = index1;
        m_last_index2  = index2;
        m_last_breakenergy  = breakenergy;

        // Compute and store value
        double eng        = energy.MeV();
        m_last_e_norm     = eng / m_last_breakenergy;
        m_last_log_e_norm = std::log(m_last_e_norm);

        if (eng < m_last_breakenergy ) {
            m_last_power      = std::pow(m_last_e_norm, m_last_index1);
        }
        else {
            m_last_power      = std::pow(m_last_e_norm, m_last_index2);
        }

    } // endif: recomputation was required

    // Return
    return;
}


/***********************************************************************//**
 * @brief Update Monte Carlo pre computation cache
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 *
 * Updates the precomputation cache for Monte Carlo simulations.
 ***************************************************************************/
void GModelSpectralBrokenPlaw::update_mc_cache(const GEnergy& emin,
                                         const GEnergy& emax) const

{
    // Case A: Index is not -1
    if (index1() != -1.0) {

        // Change in energy boundaries?
        if (emin.MeV() != m_mc_emin || emax.MeV() != m_mc_emax) {
            m_mc_emin       = emin.MeV();
            m_mc_emax       = emax.MeV();
            m_mc_exponent1   = index1() + 1.0;
            m_mc_exponent2   = index2() + 1.0;
            m_mc_pow_emin   = std::pow(m_mc_emin, m_mc_exponent1);
            m_mc_pow_ewidth = std::pow(m_mc_emax, m_mc_exponent1) - m_mc_pow_emin;
        }

    }

    // Case B: Index is -1
    else {

        // Change in energy boundaries?
        if (emin.MeV() != m_mc_emin || emax.MeV() != m_mc_emax) {
            m_mc_emin       = emin.MeV();
            m_mc_emax       = emax.MeV();
            m_mc_exponent1   = 0.0;
            m_mc_exponent2   = index2() + 1.0;
            m_mc_pow_emin   = std::log(m_mc_emin);
            m_mc_pow_ewidth = std::log(m_mc_emax) - m_mc_pow_emin;
        }

    }

    // Return
    return;
}
