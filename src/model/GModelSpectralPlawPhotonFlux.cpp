/***************************************************************************
 *   GModelSpectralPlawPhotonFlux.cpp - Spectral power law model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2016 by Juergen Knoedlseder                         *
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
 * @file GModelSpectralPlawPhotonFlux.cpp
 * @brief Flux normalized power law spectral model class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GRan.hpp"
#include "GModelSpectralPlawPhotonFlux.hpp"
#include "GModelSpectralRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpectralPlawPhotonFlux    g_spectral_plaw_phflux_seed1("PowerLaw",
                                                    "PhotonFlux",
                                                    "Index",
                                                    "LowerLimit",
                                                    "UpperLimit");
const GModelSpectralRegistry g_spectral_plaw_phflux_registry1(&g_spectral_plaw_phflux_seed1);
#if defined(G_LEGACY_XML_FORMAT)
const GModelSpectralPlawPhotonFlux    g_spectral_plaw_phflux_seed2("PowerLaw2",
                                                    "Integral",
                                                    "Index",
                                                    "LowerLimit",
                                                    "UpperLimit");
const GModelSpectralRegistry g_spectral_plaw_phflux_registry2(&g_spectral_plaw_phflux_seed2);
#endif

/* __ Method name definitions ____________________________________________ */
#define G_FLUX                "GModelSpectralPlawPhotonFlux::flux(GEnergy&, GEnergy&)"
#define G_EFLUX              "GModelSpectralPlawPhotonFlux::eflux(GEnergy&, GEnergy&)"
#define G_MC     "GModelSpectralPlawPhotonFlux::mc(GEnergy&, GEnergy&, GTime&, GRan&)"
#define G_READ                      "GModelSpectralPlawPhotonFlux::read(GXmlElement&)"
#define G_WRITE                    "GModelSpectralPlawPhotonFlux::write(GXmlElement&)"

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
 *
 * Constructs empty power law photon flux model.
 ***************************************************************************/
GModelSpectralPlawPhotonFlux::GModelSpectralPlawPhotonFlux(void) : GModelSpectral()
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
 * @param[in] photon_flux Name of photon flux parameter.
 * @param[in] index Name of index parameter.
 * @param[in] emin Name of emin parameter.
 * @param[in] emax Name of emax parameter.
 ***************************************************************************/
GModelSpectralPlawPhotonFlux::GModelSpectralPlawPhotonFlux(const std::string& type,
                                         const std::string& photon_flux,
                                         const std::string& index,
                                         const std::string& emin,
                                         const std::string& emax) :
                     GModelSpectral()
{
    // Initialise members
    init_members();

    // Set model type
    m_type = type;

    // Set parameter names
    m_photon_flux.name(photon_flux);
    m_index.name(index);
    m_emin.name(emin);
    m_emax.name(emax);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] photon_flux Photon flux (ph/cm2/s).
 * @param[in] index Power law index.
 * @param[in] emin Minimum energy.
 * @param[in] emax Maximum energy.
 *
 * Construct a spectral power law from the
 * - integral photon flux (in ph/cm2/s),
 * - spectral index,
 * - minimum energy and
 * - maximum energy.
 ***************************************************************************/
GModelSpectralPlawPhotonFlux::GModelSpectralPlawPhotonFlux(const double& photon_flux,
                                         const double&  index,
                                         const GEnergy& emin,
                                         const GEnergy& emax) :
                     GModelSpectral()
{
    // Initialise members
    init_members();

    // Set parameters
    m_photon_flux.value(photon_flux);
    m_index.value(index);
    m_emin.value(emin.MeV());
    m_emax.value(emax.MeV());

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs flux normalized power law spectral model by extracting
 * information from an XML element. See the read() method for more
 * information about the expected structure of the XML element.
 ***************************************************************************/
GModelSpectralPlawPhotonFlux::GModelSpectralPlawPhotonFlux(const GXmlElement& xml) :
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
 * @param[in] model Spectral power law model.
 ***************************************************************************/
GModelSpectralPlawPhotonFlux::GModelSpectralPlawPhotonFlux(const GModelSpectralPlawPhotonFlux& model) :
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
GModelSpectralPlawPhotonFlux::~GModelSpectralPlawPhotonFlux(void)
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
 * @param[in] model Spectral power law model.
 * @return Spectral power law model.
 ***************************************************************************/
GModelSpectralPlawPhotonFlux& GModelSpectralPlawPhotonFlux::operator=(const GModelSpectralPlawPhotonFlux& model)
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
 * @brief Clear power law model
 ***************************************************************************/
void GModelSpectralPlawPhotonFlux::clear(void)
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
 * @brief Clone power law model
 *
 * @return Pointer to deep copy of power law model
 ***************************************************************************/
GModelSpectralPlawPhotonFlux* GModelSpectralPlawPhotonFlux::clone(void) const
{
    // Clone power law model
    return new GModelSpectralPlawPhotonFlux(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @return Model value (ph/cm2/s/MeV).
 *
 * Computes
 *
 * \f[
 *    S_{\rm E}(E | t) = {\tt m\_photon\_flux}
 *    \frac{{\tt m\_index}+1}
 *         {{\tt e\_max}^{{\tt m\_index}+1} -
 *          {\tt e\_min}^{{\tt m\_index}+1}}
 *    E^{\tt m\_index}
 * \f]
 *
 * for \f${\tt m\_index} \ne -1\f$ and
 *
 * \f[
 *    S_{\rm E}(E | t) = 
 *    \frac{{\tt m\_photon\_flux}}
 *         {\log {\tt e\_max} - \log {\tt e\_min}}
 *    E^{\tt m\_index}
 * \f]
 *
 * for \f${\tt m\_index} = -1\f$, where
 * - \f${\tt e\_min}\f$ is the minimum energy of an interval,
 * - \f${\tt e\_max}\f$ is the maximum energy of an interval,
 * - \f${\tt m\_photon\_flux}\f$ is the photon flux flux between
 *   \f${\tt e\_min}\f$ and \f${\tt e\_max}\f$, and
 * - \f${\tt m\_index}\f$ is the spectral index.
 ***************************************************************************/
double GModelSpectralPlawPhotonFlux::eval(const GEnergy& srcEng,
                                 const GTime&   srcTime) const
{
    // Update precomputed values
    update(srcEng);

    // Compute function value
    double value = m_photon_flux.value() * m_norm * m_power;

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GModelSpectralPlawPhotonFlux::eval";
        std::cout << "(srcEng=" << srcEng;
        std::cout << ", srcTime=" << srcTime << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", photon_flux=" << photon_flux();
        std::cout << ", m_norm=" << m_norm;
        std::cout << ", m_power=" << m_power;
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
 * Computes
 *
 * \f[
 *    S_{\rm E}(E | t) = {\tt m\_photon\_flux}
 *    \frac{{\tt m\_index}+1}
 *         {{\tt e\_max}^{{\tt m\_index}+1} -
 *          {\tt e\_min}^{{\tt m\_index}+1}}
 *    E^{\tt m\_index}
 * \f]
 *
 * for \f${\tt m\_index} \ne -1\f$ and
 *
 * \f[
 *    S_{\rm E}(E | t) = 
 *    \frac{{\ttm\_photon\_flux}}
 *         {\log {\tt e\_max} - \log {\tt e\_min}}
 *    E^{\tt m\_index}
 * \f]
 *
 * for \f${\tt m\_index} = -1\f$, where
 * - \f${\tt e\_min}\f$ is the minimum energy of an interval,
 * - \f${\tt e\_max}\f$ is the maximum energy of an interval,
 * - \f${\tt m\_photon\_flux}\f$ is the photon flux between
 *   \f${\tt e\_min}\f$ and \f${\tt e\_max}\f$, and
 * - \f${\tt m\_index}\f$ is the spectral index.
 *
 * The method also evaluates the partial derivatives of the model with
 * respect to the parameters using
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta {\tt m\_photon\_flux}} =
 *      \frac{S_{\rm E}(E | t)}{{\tt m\_photon\_flux}}
 * \f]
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta {\tt m\_index}} =
 *      S_{\rm E}(E | t) \,
 *      \left( \frac{1}{{\tt m\_index}+1} -
 *             \frac{\log({\tt e\_max}) {\tt e\_max}^{{\tt m\_index}+1} -
 *                   \log({\tt e\_min}) {\tt e\_min}^{{\tt m\_index}+1}}
 *                        {{\tt e\_max}^{{\tt m\_index}+1} -
 *                         {\tt e\_min}^{{\tt m\_index}+1}}
 *                 + \ln(E) \right)
 * \f]
 *
 * for \f${\tt m\_index} \ne -1\f$ and
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta {\tt m\_index}} =
 *      S_{\rm E}(E | t) \, \ln(E)
 * \f]
 *
 * for \f${\tt m\_index} = -1\f$.
 *
 * No partial derivatives are supported for the energy boundaries.
 ***************************************************************************/
double GModelSpectralPlawPhotonFlux::eval_gradients(const GEnergy& srcEng,
                                           const GTime&   srcTime)
{
    // Initialise gradients
    double g_photon_flux = 0.0;
    double g_index    = 0.0;
    
    // Update precomputed values
    update(srcEng);

    // Compute function value
    double value = m_photon_flux.value() * m_norm * m_power;

    // Photon flux gradient
    if (m_photon_flux.is_free() && m_photon_flux.factor_value() > 0.0) {
         g_photon_flux = value / m_photon_flux.factor_value();
    }

    // Index gradient
    if (m_index.is_free()) {
        g_index = value * (m_g_norm + gammalib::ln10 * srcEng.log10MeV()) *
                  m_index.scale();
    }

    // Set gradients
    m_photon_flux.factor_gradient(g_photon_flux);
    m_index.factor_gradient(g_index);

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GModelSpectralPlawPhotonFlux::eval_gradients";
        std::cout << "(srcEng=" << srcEng;
        std::cout << ", srcTime=" << srcTime << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", photon_flux=" << photon_flux();
        std::cout << ", m_norm=" << m_norm;
        std::cout << ", m_power=" << m_power;
        std::cout << ", g_photon_flux=" << g_photon_flux;
        std::cout << ", g_index=" << g_index;
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
 * @param[in] emax Minimum photon energy.
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
 ***************************************************************************/
double GModelSpectralPlawPhotonFlux::flux(const GEnergy& emin,
                                 const GEnergy& emax) const
{
    // Initialise flux
    double flux = 0.0;

    // Compute only if integration range is valid
    if (emin < emax) {

        // Case A: Index is not -1
        if (index() != -1.0) {
            double gamma        = m_index.value() + 1.0;
            double pow_emin     = std::pow(emin.MeV(), gamma);
            double pow_emax     = std::pow(emax.MeV(), gamma);
            double pow_ref_emin = std::pow(this->emin().MeV(), gamma);
            double pow_ref_emax = std::pow(this->emax().MeV(), gamma);
            double factor       = (pow_emax - pow_emin) /
                                  (pow_ref_emax - pow_ref_emin);
            flux                = m_photon_flux.value() * factor;
        }

        // Case B: Index is -1
        else {
            double log_emin     = std::log(emin.MeV());
            double log_emax     = std::log(emax.MeV());
            double log_ref_emin = std::log(this->emin().MeV());
            double log_ref_emax = std::log(this->emax().MeV());
            double factor       = (log_emax - log_emin) /
                                  (log_ref_emax - log_ref_emin);
            flux                = m_photon_flux.value() * factor;
        }

    } // endif: integration range was valid

    // Return flux
    return flux;
}


/***********************************************************************//**
 * @brief Returns model energy flux between [emin, emax] (units: ph/cm2/s)
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Minimum photon energy.
 * @return Photon flux (ph/cm2/s).
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
 ***************************************************************************/
double GModelSpectralPlawPhotonFlux::eflux(const GEnergy& emin,
                                  const GEnergy& emax) const
{
    // Initialise flux
    double eflux = 0.0;

    // Compute only if integration range is valid
    if (emin < emax) {

        // Compute power law normalization
        double norm;
        if (index() != -1.0) {
            double gamma        = m_index.value() + 1.0;
            double pow_ref_emin = std::pow(this->emin().MeV(), gamma);
            double pow_ref_emax = std::pow(this->emax().MeV(), gamma);
            norm = m_photon_flux.value() * gamma / (pow_ref_emax - pow_ref_emin);
        }
        else {
            double log_ref_emin = std::log(this->emin().MeV());
            double log_ref_emax = std::log(this->emax().MeV());
            norm = m_photon_flux.value() / (log_ref_emax - log_ref_emin);
        }

        // Compute energy flux
        if (index() != -2.0) {
            double gamma    = m_index.value() + 2.0;
            double pow_emin = std::pow(emin.MeV(), gamma);
            double pow_emax = std::pow(emax.MeV(), gamma);
            eflux = norm / gamma * (pow_emax - pow_emin);
        }

        // Case B: Index is -2
        else {
            double log_emin = std::log(emin.MeV());
            double log_emax = std::log(emax.MeV());
            eflux = norm * (log_emax - log_emin);
        }

        // Convert from MeV/cm2/s to erg/cm2/s
        eflux *= gammalib::MeV2erg;

    } // endif: integration range was valid

    // Return flux
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
 * Returns Monte Carlo energy by randomly drawing from a power law.
 ***************************************************************************/
GEnergy GModelSpectralPlawPhotonFlux::mc(const GEnergy& emin,
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

    // Case A: Index is not -1
    if (index() != -1.0) {
        double exponent = index() + 1.0;
        double e_max    = std::pow(emax.MeV(), exponent);
        double e_min    = std::pow(emin.MeV(), exponent);
        double u        = ran.uniform();
        double eng      = (u > 0.0) 
                          ? std::exp(std::log(u * (e_max - e_min) + e_min) / exponent)
                          : 0.0;
        energy.MeV(eng);
    }

    // Case B: Index is -1
    else {
        double e_max = std::log(emax.MeV());
        double e_min = std::log(emin.MeV());
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
void GModelSpectralPlawPhotonFlux::read(const GXmlElement& xml)
{
    // Get parameter pointers
    const GXmlElement* phflux = gammalib::xml_get_par(G_READ, xml, m_photon_flux.name());
    const GXmlElement* index  = gammalib::xml_get_par(G_READ, xml, m_index.name());
    const GXmlElement* emin   = gammalib::xml_get_par(G_READ, xml, m_emin.name());
    const GXmlElement* emax   = gammalib::xml_get_par(G_READ, xml, m_emax.name());

    // Read parameters
    m_photon_flux.read(*phflux);
    m_index.read(*index);
    m_emin.read(*emin);
    m_emax.read(*emax);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::model_invalid_spectral
 *            Existing XML element is not of type "PowerLaw"
 *
 * Writes the spectral information into an XML element.
 ***************************************************************************/
void GModelSpectralPlawPhotonFlux::write(GXmlElement& xml) const
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
    GXmlElement* phflux = gammalib::xml_need_par(G_WRITE, xml, m_photon_flux.name());
    GXmlElement* index  = gammalib::xml_need_par(G_WRITE, xml, m_index.name());
    GXmlElement* emin   = gammalib::xml_need_par(G_WRITE, xml, m_emin.name());
    GXmlElement* emax   = gammalib::xml_need_par(G_WRITE, xml, m_emax.name());

    // Write parameters
    m_photon_flux.write(*phflux);
    m_index.write(*index);
    m_emin.write(*emin);
    m_emax.write(*emax);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print power law information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing power law information.
 ***************************************************************************/
std::string GModelSpectralPlawPhotonFlux::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpectralPlawPhotonFlux ===");

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
void GModelSpectralPlawPhotonFlux::init_members(void)
{
    // Initialise model type
    m_type = "PowerLaw";

    // Initialise photon flux
    m_photon_flux.clear();
    m_photon_flux.name("PhotonFlux");
    m_photon_flux.unit("ph/cm2/s");
    m_photon_flux.scale(1.0);
    m_photon_flux.value(1.0);       // default: 1.0
    m_photon_flux.range(0.0, 10.0); // range:   [0,10]
    m_photon_flux.free();
    m_photon_flux.gradient(0.0);
    m_photon_flux.has_grad(true);

    // Initialise powerlaw index
    m_index.clear();
    m_index.name("Index");
    m_index.scale(1.0);
    m_index.value(-2.0);        // default: -2.0
    m_index.range(-10.0,+10.0); // range:   [-10,+10]
    m_index.free();
    m_index.gradient(0.0);
    m_index.has_grad(true);

    // Initialise lower limit
    m_emin.clear();
    m_emin.name("LowerLimit");
    m_emin.unit("MeV");
    m_emin.scale(1.0);
    m_emin.value(100.0);         // default: 100
    m_emin.range(0.001, 1.0e15); // range:   [0.001, 1e15]
    m_emin.fix();
    m_emin.gradient(0.0);
    m_emin.has_grad(false);

    // Initialise upper limit
    m_emax.clear();
    m_emax.name("UpperLimit");
    m_emax.unit("MeV");
    m_emax.scale(1.0);
    m_emax.value(500000.0);      // default: 500000
    m_emin.range(0.001, 1.0e15); // range:   [0.001, 1e15]
    m_emax.fix();
    m_emax.gradient(0.0);
    m_emax.has_grad(false);

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_photon_flux);
    m_pars.push_back(&m_index);
    m_pars.push_back(&m_emin);
    m_pars.push_back(&m_emax);

    // Initialise last parameters (for fast computation)
    m_log_emin        = 0.0;
    m_log_emax        = 0.0;
    m_pow_emin        = 0.0;
    m_pow_emax        = 0.0;
    m_norm            = 0.0;
    m_power           = 0.0;
    m_last_phflux     = 0.0;
    m_last_index      = 1000.0;
    m_last_emin.MeV(0.0);
    m_last_emax.MeV(0.0);
    m_last_energy.MeV(0.0);
    m_last_value      = 0.0;
    m_last_g_phflux   = 0.0;
    m_last_g_index    = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Spectral power law model.
 ***************************************************************************/
void GModelSpectralPlawPhotonFlux::copy_members(const GModelSpectralPlawPhotonFlux& model)
{
    // Copy members
    m_type        = model.m_type;
    m_photon_flux = model.m_photon_flux;
    m_index       = model.m_index;
    m_emin        = model.m_emin;
    m_emax        = model.m_emax;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_photon_flux);
    m_pars.push_back(&m_index);
    m_pars.push_back(&m_emin);
    m_pars.push_back(&m_emax);

    // Copy bookkeeping information
    m_log_emin        = model.m_log_emin;
    m_log_emax        = model.m_log_emax;
    m_pow_emin        = model.m_pow_emin;
    m_pow_emax        = model.m_pow_emax;
    m_norm            = model.m_norm;
    m_power           = model.m_power;
    m_last_phflux   = model.m_last_phflux;
    m_last_index      = model.m_last_index;
    m_last_emin       = model.m_last_emin;
    m_last_emax       = model.m_last_emax;
    m_last_energy     = model.m_last_energy;
    m_last_value      = model.m_last_value;
    m_last_g_phflux   = model.m_last_g_phflux;
    m_last_g_index    = model.m_last_g_index;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralPlawPhotonFlux::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update precomputed values
 *
 * @param[in] srcEng Source energy
 ***************************************************************************/
void GModelSpectralPlawPhotonFlux::update(const GEnergy& srcEng) const
{
    // Compute index+1
    double gamma = index() + 1.0;

    // Change in spectral index?
    if (index() != m_last_index) {

        // Save actual spectral index
        m_last_index = index();

        // Change in energy boundaries?
        if (emin() != m_last_emin || emax() != m_last_emax) {
            m_log_emin  = std::log(emin().MeV());
            m_last_emin = emin();
            m_log_emax  = std::log(emax().MeV());
            m_last_emax = emax();
        }

        // Compute normalization factors
        if (gamma != 0.0) {
            m_pow_emin  = std::pow(emin().MeV(), gamma);
            m_pow_emax  = std::pow(emax().MeV(), gamma);
            double d    = m_pow_emax - m_pow_emin;
            m_norm      = gamma / d;
            m_g_norm    = 1.0/gamma - 
                          (m_pow_emax*m_log_emax - m_pow_emin*m_log_emin)/d;
        }
        else {
            m_norm   = 1.0 / (m_log_emax - m_log_emin);
            m_g_norm = 0.0;
        }

        // Update power law factor
        m_power       = std::pow(srcEng.MeV(), index());
        m_last_energy = srcEng;

    } // endif: change in spectral index

    // ... no change in spectral index
    else {

        // Change in energy boundaries?
        if (emin() != m_last_emin || emax() != m_last_emax) {

            // Update energy boundaries
            m_log_emin  = std::log(emin().MeV());
            m_last_emin = emin();
            m_log_emax  = std::log(emax().MeV());
            m_last_emax = emax();
        
            // Compute power law normalization
            if (gamma != 0.0) {
                m_pow_emin  = std::pow(emin().MeV(), gamma);
                m_pow_emax  = std::pow(emax().MeV(), gamma);
                double d    = m_pow_emax - m_pow_emin;
                m_norm      = gamma / d;
                m_g_norm    = 1.0/gamma - 
                              (m_pow_emax*m_log_emax - m_pow_emin*m_log_emin)/d;
            }
            else {
                m_norm = 1.0 / (m_log_emax - m_log_emin);
                m_g_norm = 0.0;
            }

        } // endif: change in energy boundaries

        // Change in energy?
        if (srcEng != m_last_energy) {
            m_power       = std::pow(srcEng.MeV(), index());
            m_last_energy = srcEng;
        }

    } // endelse: no change in spectral index

    // Return
    return;
}
