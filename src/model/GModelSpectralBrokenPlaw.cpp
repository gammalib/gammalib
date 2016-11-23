/***************************************************************************
 *     GModelSpectralBrokenPlaw.cpp - Broken power law spectrum class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2016 by Anneli Schulz                               *
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
 * @brief Broken power law spectrum class implementation
 * @author Anneli Schulz
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GRan.hpp"
#include "GModelSpectralBrokenPlaw.hpp"
#include "GModelSpectralRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpectralBrokenPlaw g_spectral_blaw_seed1("BrokenPowerLaw",
                                                     "Prefactor",
                                                     "Index1",
                                                     "BreakEnergy",
                                                     "Index2");
const GModelSpectralRegistry   g_spectral_blaw_registry1(&g_spectral_blaw_seed1);
#if defined(G_LEGACY_XML_FORMAT)
const GModelSpectralBrokenPlaw g_spectral_blaw_seed2("BrokenPowerLaw",
                                                     "Prefactor",
                                                     "Index1",
                                                     "BreakValue",
                                                     "Index2");
const GModelSpectralRegistry   g_spectral_blaw_registry2(&g_spectral_blaw_seed2);
#endif

/* __ Method name definitions ____________________________________________ */
#define G_MC       "GModelSpectralBrokenPlaw::mc(GEnergy&, GEnergy&, GTime&,"\
                                                                    " GRan&)"
#define G_READ                 "GModelSpectralBrokenPlaw::read(GXmlElement&)"
#define G_WRITE               "GModelSpectralBrokenPlaw::write(GXmlElement&)"

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
 * @brief Model type and parameter name constructor
 *
 * @param[in] type Model type.
 * @param[in] prefactor Name of prefactor parameter.
 * @param[in] index1 Name of index1 parameter.
 * @param[in] breakenergy Name of breakenergy parameter.
 * @param[in] index2 Name of index2 parameter.
 ***************************************************************************/
GModelSpectralBrokenPlaw::GModelSpectralBrokenPlaw(const std::string& type,
                                                   const std::string& prefactor,
                                                   const std::string& index1,
                                                   const std::string& breakenergy,
                                                   const std::string& index2) :
                          GModelSpectral()
{
    // Initialise members
    init_members();

    // Set model type
    m_type = type;

    // Set parameter names
    m_norm.name(prefactor);
    m_index1.name(index1);
    m_breakenergy.name(breakenergy);
    m_index2.name(index2);

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
 *
 *     power law @p prefactor (ph/cm2/s/MeV)
 *     spectral @p index1
 *     Energy @p breakenergy of spectral break
 *     spectral @p index2
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
GModelSpectralBrokenPlaw::GModelSpectralBrokenPlaw(const GXmlElement& xml) :
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
 * @param[in] model Broken power law model.
 ***************************************************************************/
GModelSpectralBrokenPlaw::GModelSpectralBrokenPlaw(const GModelSpectralBrokenPlaw& model) :
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
 * @param[in] model Broken power law model.
 * @return Broken power law model.
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
 * @brief Clear broken power law model
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
 * @brief Clone broken power law model
 *
 * @return Pointer to deep copy of broken power law model.
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
 * @param[in] gradients Compute gradients?
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
 * If the @p gradients flag is true the method will also compute the
 * partial derivatives of the model with respect to the parameters using
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
double GModelSpectralBrokenPlaw::eval(const GEnergy& srcEng,
                                      const GTime&   srcTime,
                                      const bool&    gradients) const
{
    // Update the evaluation cache
    update_eval_cache(srcEng);

    // Compute function value
    double value = m_norm.value() * m_last_power;

    // Optionally compute gradients
    if (gradients) {

        // Compute normalisation gradient
        double g_norm  = (m_norm.is_free()) ? m_norm.scale() * m_last_power : 0.0;
        m_norm.factor_gradient(g_norm);

        // Compute index and break value gradients
        if (srcEng.MeV() < m_breakenergy.value()) {
            double g_index = (m_index1.is_free())
                             ? value * m_index1.scale() * m_last_log_e_norm
                             : 0.0;
            double g_break = (m_breakenergy.is_free())
                             ? -value * m_last_index1 / m_breakenergy.factor_value()
                             : 0.0;
            m_index1.factor_gradient(g_index);
            m_index2.factor_gradient(0.0);
            m_breakenergy.factor_gradient(g_break);
        }
        else {
            double g_index = (m_index2.is_free())
                             ? value * m_index2.scale() * m_last_log_e_norm
                             : 0.0;
            double g_break = (m_breakenergy.is_free())
                             ? -value * m_last_index2 / m_breakenergy.factor_value()
                             : 0.0;
            m_index1.factor_gradient(0.0);
            m_index2.factor_gradient(g_index);
            m_breakenergy.factor_gradient(g_break);
        }

    } // endif: gradient computation was requested

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
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
        
        // First case: emin < breakenergy < emax
        if (emin.MeV() < m_breakenergy.value() && m_breakenergy.value() < emax.MeV()) {

            // Compute photon flux
            flux = m_norm.value() *
                   (gammalib::plaw_photon_flux(emin.MeV(),
                                               m_breakenergy.value(),
                                               m_breakenergy.value(),
                                               m_index1.value()) +
                    gammalib::plaw_photon_flux(m_breakenergy.value(),
                                               emax.MeV(),
                                               m_breakenergy.value(),
                                               m_index2.value()));
        }
        
        // Second case: breakenergy > emax:
        else if (m_breakenergy.value() > emax.MeV()) {
            flux = m_norm.value() *
                   gammalib::plaw_photon_flux(emin.MeV(),
                                              emax.MeV(),
                                              m_breakenergy.value(),
                                              m_index1.value());
        }
        
        // Third case breakenergy < emin:
        else if (m_breakenergy.value() < emin.MeV()) {
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

        // First case: emin < breakenergy < emax
        if (emin.MeV() < m_breakenergy.value() && m_breakenergy.value() < emax.MeV()) {
            eflux = m_norm.value() *
                   (gammalib::plaw_energy_flux(emin.MeV(),
                                               m_breakenergy.value(),
                                               m_breakenergy.value(),
                                               m_index1.value()) +
                    gammalib::plaw_energy_flux(m_breakenergy.value(),
                                               emax.MeV(),
                                               m_breakenergy.value(),
                                               m_index2.value()));
        }

        // Second case: breakenergy > emax:
        else if (m_breakenergy.value() > emax.MeV()) {
            eflux = m_norm.value() *
                   gammalib::plaw_energy_flux(emin.MeV(),
                                              emax.MeV(),
                                              m_breakenergy.value(),
                                              m_index1.value());
        }

        // Third case breakenergy < emin:
        else if (m_breakenergy.value() < emin.MeV()) {
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

    // Allocate energy
    GEnergy energy;

    // Update cache
    update_mc_cache(emin, emax);

    // Determine in which bin we reside
    int inx = 0;
    if (m_mc_cum.size() > 1) {
        double u = ran.uniform();
        for (inx = m_mc_cum.size()-1; inx > 0; --inx) {
            if (m_mc_cum[inx-1] <= u)
                break;
        }
    }

    // Get random energy for specific bin
    if (m_mc_exp[inx] != 0.0) {
        double e_min = m_mc_min[inx];
        double e_max = m_mc_max[inx];
        double u     = ran.uniform();
        double eng   = (u > 0.0)
                        ? std::exp(std::log(u * (e_max - e_min) + e_min) / m_mc_exp[inx])
                        : 0.0;
        energy.MeV(eng);
    }
    else {
        double e_min = m_mc_min[inx];
        double e_max = m_mc_max[inx];
        double u     = ran.uniform();
        double eng   = std::exp(u * (e_max - e_min) + e_min);
        energy.MeV(eng);
    }

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
void GModelSpectralBrokenPlaw::read(const GXmlElement& xml)
{
    // Get remaining XML parameters
    const GXmlElement* prefactor   = gammalib::xml_get_par(G_READ, xml, m_norm.name());
    const GXmlElement* index1      = gammalib::xml_get_par(G_READ, xml, m_index1.name());
    const GXmlElement* breakenergy = gammalib::xml_get_par(G_READ, xml, m_breakenergy.name());
    const GXmlElement* index2      = gammalib::xml_get_par(G_READ, xml, m_index2.name());

    // Read parameters
    m_norm.read(*prefactor);
    m_index1.read(*index1);
    m_breakenergy.read(*breakenergy);
    m_index2.read(*index2);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::model_invalid_spectral
 *            Existing XML element is not of the expected type.
 *
 * Writes the spectral information into an XML element.
 ***************************************************************************/
void GModelSpectralBrokenPlaw::write(GXmlElement& xml) const
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

    // Get XML parameters
    GXmlElement* norm        = gammalib::xml_need_par(G_WRITE, xml, m_norm.name());
    GXmlElement* index1      = gammalib::xml_need_par(G_WRITE, xml, m_index1.name());
    GXmlElement* index2      = gammalib::xml_need_par(G_WRITE, xml, m_index2.name());
    GXmlElement* breakenergy = gammalib::xml_need_par(G_WRITE, xml, m_breakenergy.name());

    // Write parameters
    m_norm.write(*norm);
    m_index1.write(*index1);
    m_index2.write(*index2);
    m_breakenergy.write(*breakenergy);

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
    // Initialise model type
    m_type = "BrokenPowerLaw";

    // Initialise powerlaw normalisation
    m_norm.clear();
    m_norm.name("Prefactor");
    m_norm.unit("ph/cm2/s/MeV");
    m_norm.scale(1.0);
    m_norm.value(1.0);          // default: 1.0
    m_norm.min(0.0);            // min:     0.0
    m_norm.free();
    m_norm.gradient(0.0);
    m_norm.has_grad(true);

    // Initialise powerlaw index1
    m_index1.clear();
    m_index1.name("Index1");
    m_index1.scale(1.0);
    m_index1.value(-2.0);        // default: -2.0
    m_index1.range(-10.0,+10.0); // range:   [-10,+10]
    m_index1.free();
    m_index1.gradient(0.0);
    m_index1.has_grad(true);

    // Initialise powerlaw index2
    m_index2.clear();
    m_index2.name("Index2");
    m_index2.scale(1.0);
    m_index2.value(-2.0);        // default: -2.0
    m_index2.range(-10.0,+10.0); // range:   [-10,+10]
    m_index2.free();
    m_index2.gradient(0.0);
    m_index2.has_grad(true);

    // Initialise break energy
    m_breakenergy.clear();
    m_breakenergy.name("BreakEnergy");
    m_breakenergy.unit("MeV");
    m_breakenergy.scale(1.0);
    m_breakenergy.value(100.0);  // default: 100
    m_breakenergy.fix();
    m_breakenergy.gradient(0.0);
    m_breakenergy.has_grad(true);

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
    m_last_breakenergy = 1.0e30;
    m_last_e_norm      = 0.0;
    m_last_log_e_norm  = 0.0;
    m_last_power       = 0.0;

    // Initialise MC cache
    m_mc_exponent1  = 0.0;
    m_mc_exponent2  = 0.0;
    m_mc_pow_emin   = 0.0;
    m_mc_pow_ewidth = 0.0;

    // Initialise MC cache
    m_mc_emin = 0.0;
    m_mc_emax = 0.0;


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
    m_type        = model.m_type;
    m_norm        = model.m_norm;
    m_index1      = model.m_index1;
    m_breakenergy = model.m_breakenergy;
    m_index2      = model.m_index2;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);
    m_pars.push_back(&m_index1);
    m_pars.push_back(&m_breakenergy);
    m_pars.push_back(&m_index2);

    // Copy eval cache
    m_last_energy      = model.m_last_energy;
    m_last_index1      = model.m_last_index1;
    m_last_breakenergy = model.m_last_breakenergy;
    m_last_index2      = model.m_last_index2;
    m_last_e_norm      = model.m_last_e_norm;
    m_last_log_e_norm  = model.m_last_log_e_norm;
    m_last_power       = model.m_last_power;

    // Copy MC cache
    m_mc_emin       = model.m_mc_emin;
    m_mc_emax       = model.m_mc_emax;
    m_mc_exponent1  = model.m_mc_exponent1;
    m_mc_exponent2  = model.m_mc_exponent2;
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
 * Updates the precomputation cache for eval() the method.
 ***************************************************************************/
void GModelSpectralBrokenPlaw::update_eval_cache(const GEnergy& energy) const
{
    // Get parameter values (takes 2 multiplications which are difficult
    // to avoid)
    double index1      = m_index1.value();
    double breakenergy = m_breakenergy.value();
    double index2      = m_index2.value();

    // If the energy or one of the parameters index1, index2 or breakenergy
    // energy has changed then recompute the cache
    if ((m_last_energy       != energy) ||
        (m_last_index1       != index1)  ||
        (m_last_index2       != index2)  ||
        (m_last_breakenergy  != breakenergy)) {

        // Store actual energy and parameter values
        m_last_energy       = energy;
        m_last_index1       = index1;
        m_last_index2       = index2;
        m_last_breakenergy  = breakenergy;

        // Compute and store value
        double eng        = energy.MeV();
        m_last_e_norm     = eng / m_last_breakenergy;
        m_last_log_e_norm = std::log(m_last_e_norm);
        if (eng < m_last_breakenergy ) {
            m_last_power = std::pow(m_last_e_norm, m_last_index1);
        }
        else {
            m_last_power = std::pow(m_last_e_norm, m_last_index2);
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
    // Check if we need to update the cache
    if (emin.MeV() != m_mc_emin || emax.MeV() != m_mc_emax) {

        // Store new energy interval
        m_mc_emin = emin.MeV();
        m_mc_emax = emax.MeV();

        // Initialise cache
        m_mc_cum.clear();
        m_mc_min.clear();
        m_mc_max.clear();
        m_mc_exp.clear();

        // Get energy range in MeV
        double e_min = emin.MeV();
        double e_max = emax.MeV();

        // Continue only if e_max > e_min
        if (e_max > e_min) {

            // Allocate flux
            double flux;

            // Determine left node index for minimum and maximum energy
            int inx_emin = (e_min < m_breakenergy.value() ) ? 0 : 1;
            int inx_emax = (e_max < m_breakenergy.value() ) ? 0 : 1;

            // If both energies are within the same node then just
            // add this one node on the stack
            if (inx_emin == inx_emax) {
                double exp_valid = (e_min < m_breakenergy.value())
                                   ? m_index1.value() : m_index2.value();
                flux = m_norm.value() *
                       gammalib::plaw_photon_flux(e_min,
                                                  e_max,
                                                  m_breakenergy.value(),
                                                  exp_valid);
                m_mc_cum.push_back(flux);
                m_mc_min.push_back(e_min);
                m_mc_max.push_back(e_max);
                m_mc_exp.push_back(exp_valid);
            }

            // ... otherwise integrate over both nodes
            else {
                // just enter the values for first pl: bin [0]
                flux = m_norm.value() *
                       gammalib::plaw_photon_flux(e_min,
                                                  m_breakenergy.value(),
                                                  m_breakenergy.value(),
                                                  m_index1.value());
                m_mc_cum.push_back(flux);
                m_mc_exp.push_back(m_index1.value());
                m_mc_min.push_back(e_min);
                m_mc_max.push_back(m_breakenergy.value());

                // and for
                flux = m_norm.value() *
                       gammalib::plaw_photon_flux(m_breakenergy.value(),
                                                  e_max,
                                                  m_breakenergy.value(),
                                                  m_index2.value());
                m_mc_cum.push_back(flux);
                m_mc_exp.push_back(m_index2.value());
                m_mc_max.push_back(e_max);
                m_mc_min.push_back(m_breakenergy.value());
            } // endelse: emin and emax not between same nodes

            // Build cumulative distribution
            for (int i = 1; i < m_mc_cum.size(); ++i) {
                m_mc_cum[i] += m_mc_cum[i-1];
            }
            double norm = m_mc_cum[m_mc_cum.size()-1];
            for (int i = 0; i < m_mc_cum.size(); ++i) {
                m_mc_cum[i] /= norm;
            }

            // Set MC values
            for (int i = 0; i < m_mc_cum.size(); ++i) {

                // Compute exponent
                double exponent = m_mc_exp[i] + 1.0;

                // Exponent dependend computation
                if (std::abs(exponent) > 1.0e-11) {
                    m_mc_exp[i] = exponent;
                    m_mc_min[i] = std::pow(m_mc_min[i], exponent);
                    m_mc_max[i] = std::pow(m_mc_max[i], exponent);
                }
                else {
                    m_mc_exp[i] = 0.0;
                    m_mc_min[i] = std::log(m_mc_min[i]);
                    m_mc_max[i] = std::log(m_mc_max[i]);
                }

            } // endfor: set MC values

        } // endif: e_max > e_min

    } // endif: Update was required

    // Return
    return;
}
