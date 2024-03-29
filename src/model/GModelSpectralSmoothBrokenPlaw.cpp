/***************************************************************************
 *                 GModelSpectralSmoothBrokenPlaw.cpp                      *
 *               Smoothly broken power law spectrum class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2021 by Joshua Cardenzana                           *
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
 * @file GModelSpectralSmoothBrokenPlaw.cpp
 * @brief Smoothly broken power law spectrum class implementation
 * @author Joshua Cardenzana
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GException.hpp"
#include "GIntegral.hpp"
#include "GTools.hpp"
#include "GRan.hpp"
#include "GModelSpectralSmoothBrokenPlaw.hpp"
#include "GModelSpectralRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpectralSmoothBrokenPlaw g_spectral_sblaw_seed1("SmoothBrokenPowerLaw",
                                                            "Prefactor",
                                                            "Index1",
                                                            "PivotEnergy",
                                                            "Index2",
                                                            "BreakEnergy",
                                                            "BreakSmoothness");
const GModelSpectralRegistry         g_spectral_sblaw_registry1(&g_spectral_sblaw_seed1);
const GModelSpectralSmoothBrokenPlaw g_spectral_sblaw_seed2("SmoothBrokenPowerLaw",
                                                            "Prefactor",
                                                            "Index1",
                                                            "Scale",
                                                            "Index2",
                                                            "BreakValue",
                                                            "Beta");
const GModelSpectralRegistry         g_spectral_sblaw_registry2(&g_spectral_sblaw_seed2);

/* __ Method name definitions ____________________________________________ */
#define G_MC "GModelSpectralSmoothBrokenPlaw::mc(GEnergy&, GEnergy&, GTime&,"\
                                                                    " GRan&)"
#define G_READ           "GModelSpectralSmoothBrokenPlaw::read(GXmlElement&)"
#define G_WRITE         "GModelSpectralSmoothBrokenPlaw::write(GXmlElement&)"

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
GModelSpectralSmoothBrokenPlaw::GModelSpectralSmoothBrokenPlaw(void) :
                                GModelSpectral()
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
 * @param[in] pivot Name of pivot parameter.
 * @param[in] index2 Name of index2 parameter.
 * @param[in] breakenergy Name of breakenergy parameter.
 * @param[in] beta Name of beta parameter.
 ***************************************************************************/
GModelSpectralSmoothBrokenPlaw::GModelSpectralSmoothBrokenPlaw(
                                        const std::string& type,
                                        const std::string& prefactor,
                                        const std::string& index1,
                                        const std::string& pivot,
                                        const std::string& index2,
                                        const std::string& breakenergy,
                                        const std::string& beta) :
                                GModelSpectral()
{
    // Initialise members
    init_members();
    
    // Set model type
    m_type = type;
    
    // Set parameter names
    m_norm.name(prefactor);
    m_index1.name(index1);
    m_pivot.name(pivot);
    m_index2.name(index2);
    m_breakenergy.name(breakenergy);
    m_beta.name(beta);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Parameter constructor
 *
 * @param[in] prefactor Smoothly broken power law pre factor (ph/cm2/s/MeV).
 * @param[in] index1 Smoothly broken power law index1.
 * @param[in] pivot Smoothly broken power law pivot energy
 * @param[in] index2 Smoothly broken power law index1.
 * @param[in] breakenergy Break energy.
 * @param[in] beta Break smoothness parameter
 *
 * Constructs a smoothly broken power law using the model parameters
 * power law @p prefactor (ph/cm2/s/MeV),
 * spectral @p index1,
 * @p pivot energy,
 * spectral @p index2,
 * @p breakenergy of spectral break, and
 * smoothness parameter @p beta.
 ***************************************************************************/
GModelSpectralSmoothBrokenPlaw::GModelSpectralSmoothBrokenPlaw(
                                        const double&  prefactor,
                                        const double&  index1,
                                        const GEnergy& pivot,
                                        const double&  index2,
                                        const GEnergy& breakenergy,
                                        const double&  beta) :
                                GModelSpectral()
{
    // Initialise members
    init_members();
    
    // Set parameters
    m_norm.value(prefactor);
    m_index1.value(index1);
    m_pivot.value(pivot.MeV());              // Internally stored in MeV
    m_index2.value(index2);
    m_breakenergy.value(breakenergy.MeV());  // Internally stored in MeV
    m_beta.value(beta);
    
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
 * Constructs a smoothly broken power law spectral model by extracting
 * information from an XML element. See the read() method for more
 * information about the expected structure of the XML element.
 ***************************************************************************/
GModelSpectralSmoothBrokenPlaw::GModelSpectralSmoothBrokenPlaw(const GXmlElement& xml) :
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
 * @param[in] model Smoothly broken power law model.
 ***************************************************************************/
GModelSpectralSmoothBrokenPlaw::GModelSpectralSmoothBrokenPlaw(
                                const GModelSpectralSmoothBrokenPlaw& model) :
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
GModelSpectralSmoothBrokenPlaw::~GModelSpectralSmoothBrokenPlaw(void)
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
 * @param[in] model Smoothly broken power law model.
 * @return Smoothly broken power law model.
 ***************************************************************************/
GModelSpectralSmoothBrokenPlaw& GModelSpectralSmoothBrokenPlaw::operator=(
                                const GModelSpectralSmoothBrokenPlaw& model)
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
 * @brief Clear smoothly broken power law model
 ***************************************************************************/
void GModelSpectralSmoothBrokenPlaw::clear(void)
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
 * @brief Clone smoothly broken power law model
 *
 * @return Pointer to deep copy of smoothly broken power law model.
 ***************************************************************************/
GModelSpectralSmoothBrokenPlaw* GModelSpectralSmoothBrokenPlaw::clone(void) const
{
    // Clone spectral broken power law model
    return new GModelSpectralSmoothBrokenPlaw(*this);
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
 *    S_{\rm E}(E | t) = k_0 \left( \frac{E}{E_0} \right)^{\gamma_1}
 *     \left[ 1 + \left( \frac{E}{E_b} \right)^{\frac{\gamma_1 - \gamma_2}{\beta}} 
 *     \right]^{-\beta}
 * \f]
 *
 * where:
 * - \f$k_0\f$ is the normalization or prefactor,
 * - \f$\gamma_1\f$ is the spectral index before the break,
 * - \f$\gamma_2\f$ is the spectral index after the break,
 * - \f$E_0\f$ is the pivot energy,
 * - \f$E_b\f$ is the break energy,
 * - \f$\beta\f$ is the break smoothness.
 *
 * If the @p gradients flag is true the method will also compute the
 * partial derivatives of the model with respect to the parameters using
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta k_0} =
 *      \frac{S_{\rm E}(E | t)}{k_0}
 * \f]
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta \gamma_1} =
 *      S_{\rm E}(E | t) \left[ \ln\left( \frac{E}{E_0} \right) -
 *      \frac{\left(\frac{E}{E_b}\right)^{\frac{\gamma_1 - \gamma_2}{\beta}} 
 *            \ln\left( \frac{E}{E_b} \right)}
 *           {\left(\frac{E}{E_b}\right)^{\frac{\gamma_1 - \gamma_2}{\beta}} + 1}\right]
 * \f]
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta \gamma_2} =
 *      S_{\rm E}(E | t) 
 *      \frac{\left(\frac{E}{E_b}\right)^{\frac{\gamma_1 - \gamma_2}{\beta}}
 *            \ln\left( \frac{E}{E_b} \right)}
 *           {\left(\frac{E}{E_b}\right)^{\frac{\gamma_1 - \gamma_2}{\beta}} + 1}
 * \f]
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta E_0} = 
 *      S_{\rm E}(E | t) \frac{\gamma_1}{E_0}
 * \f]
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta E_b} =
 *      S_{\rm E}(E | t) \frac{\left(\gamma_1 - \gamma_2 \right)
 *          \left( \frac{E}{E_b} \right)^{\frac{\gamma_1 - \gamma_2}{\beta}}}
        {E_b \left(1 + \left(\frac{E}{E_b}\right)^{\frac{\gamma_1 - \gamma_2}{\beta}} \right)}
 * \f]
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta \beta} = S_{\rm E}(E | t)
 *      \left[
 *         \frac{ \left( \frac{E}{E_b} \right)^{\frac{\gamma_1 - \gamma_2}{\beta}}
 *               \ln \left( \left( \frac{E}{E_b}\right)^{ \frac{\gamma_1 - \gamma_2}{\beta}} \right)}
 *              {1 + \left(\frac{E}{E_b}\right)^{\frac{\gamma_1 - \gamma_2}{\beta}}}
 *         - \ln \left( 1 + \left(\frac{E}{E_b}\right)^{\frac{\gamma_1 - \gamma_2}{\beta}} \right)
 *      \right]
 * \f]
 *
 * @todo The method expects that energy!=0. Otherwise Inf or NaN may result.
 ***************************************************************************/
double GModelSpectralSmoothBrokenPlaw::eval(const GEnergy& srcEng,
                                            const GTime&   srcTime,
                                            const bool&    gradients) const
{
    // Update the evaluation cache
    update_eval_cache(srcEng);
    
    // Compute function value
    double value = m_norm.value() * m_last_epivot_pow *
                   std::pow(1.0 + m_last_ebreak_pow, -m_last_beta);
    
    // Optionally compute gradients
    if (gradients) {
        
        // Compute normalisation gradient
        double g_norm  = (m_norm.is_free())
                ? m_norm.scale() * m_last_epivot_pow *
                    std::pow(1.0 + m_last_ebreak_pow, -m_last_beta)
                : 0.0;
        
        // Compute index1 and index2 value gradients
        double g_index1 = (m_index1.is_free())
                ? value * m_index1.scale() * (m_last_log_epivot_norm -
                  ((m_last_ebreak_pow * m_last_log_ebreak_norm) /
                   (m_last_ebreak_pow + 1.0)))
                : 0.0;
        double g_index2 = (m_index2.is_free())
                ? value * m_index2.scale() * m_last_log_ebreak_norm *
                  m_last_ebreak_pow / (1.0 + m_last_ebreak_pow)
                : 0.0;

        // Compute pivot and break energy value gradients
        double g_pivot = (m_pivot.is_free() && m_pivot.factor_value() != 0.0)
                ? -value * m_last_index1 / m_pivot.factor_value()
                : 0.0;
        double g_break = (m_breakenergy.is_free())
                ? value * (m_last_index1-m_last_index2) * m_last_ebreak_pow /
                  ((1.0+m_last_ebreak_pow) * m_breakenergy.factor_value())
                : 0.0;

        // Compute beta gradient
        double g_beta  = (m_beta.is_free())
                ? value * m_beta.scale() *
                  ((std::log(m_last_ebreak_pow) * m_last_ebreak_pow) /
                  (1.0+m_last_ebreak_pow) - std::log(1.0+m_last_ebreak_pow))
                : 0.0;
        
        // Store the gradient values
        m_norm.factor_gradient(g_norm);
        m_index1.factor_gradient(g_index1);
        m_index2.factor_gradient(g_index2);
        m_pivot.factor_gradient(g_pivot);
        m_breakenergy.factor_gradient(g_break);
        m_beta.factor_gradient(g_beta);
        
    } // endif: gradient computation was requested
    
    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GModelSpectralSmoothBrokenPlaw::eval";
        std::cout << "(srcEng=" << srcEng;
        std::cout << ", srcTime=" << srcTime << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", m_norm=" << m_norm.value();
        std::cout << ", m_index1=" << m_index1.value();
        std::cout << ", m_index2=" << m_index2.value();
        std::cout << ", m_pivot=" << m_pivot.value();
        std::cout << ", m_breakenergy=" << m_breakenergy.value();
        std::cout << ", m_beta=" << m_beta.value();
        std::cout << ", m_last_epivot_norm=" << m_last_epivot_norm;
        std::cout << ", m_last_ebreak_norm=" << m_last_ebreak_norm;
        std::cout << ", m_last_log_epivot_norm=" << m_last_log_epivot_norm;
        std::cout << ", m_last_log_ebreak_norm=" << m_last_log_ebreak_norm;
        std::cout << ", m_last_epivot_pow=" << m_last_epivot_pow;
        std::cout << ", m_last_ebreak_pow=" << m_last_ebreak_pow;
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
 * The integration is done numerically.
 ***************************************************************************/
double GModelSpectralSmoothBrokenPlaw::flux(const GEnergy& emin,
                                            const GEnergy& emax) const
{
    // Initialise flux
    double flux = 0.0;
    
    // Compute only if integration range is valid
    if (emin < emax) {
        
        // Initialise function to integrate
        flux_kern kernel(prefactor(), index1(), pivot(),
                         index2(), breakenergy(), beta());
        
        // Initialise integral class with function
        GIntegral integral(&kernel);
        
        // Set integration precision
        integral.eps(1.0e-8);
        
        // Calculate integral between emin and emax
        flux = integral.romberg(emin.MeV(), emax.MeV());
        
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
 * The integration is done numerically.
 ***************************************************************************/
double GModelSpectralSmoothBrokenPlaw::eflux(const GEnergy& emin,
                                             const GEnergy& emax) const
{
    // Initialise flux
    double eflux = 0.0;
    
    // Compute only if integration range is valid
    if (emin < emax) {
        
        // Initialise function to integrate
        eflux_kern kernel(prefactor(), index1(), pivot(),
                          index2(), breakenergy(), beta());
        
        // Initialise integral class with function
        GIntegral integral(&kernel);
        
        // Set integration precision
        integral.eps(1.0e-8);
        
        // Calculate integral between emin and emax
        eflux = integral.romberg(emin.MeV(), emax.MeV());
        
        // Convert from MeV/cm2/s to erg/cm2/s
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
 * Returns Monte Carlo energy by randomly drawing from a smoothly broken
 * power law.
 ***************************************************************************/
GEnergy GModelSpectralSmoothBrokenPlaw::mc(const GEnergy& emin,
                                           const GEnergy& emax,
                                           const GTime&   time,
                                           GRan&          ran) const
{
    // Check energy interval
    gammalib::check_energy_interval(G_MC, emin, emax);

    // Allocate energy
    GEnergy energy;
    
    // Update Monte Carlo cache
    update_mc_cache();
    
    // Initialse acceptance fraction
    double acceptance_fraction(0.0);
    
    // Use rejection method to draw a random energy. We first draw
    // analytically from a broken power law, and then compare the power law
    // at the drawn energy to the curved function. This gives an acceptance
    // fraction, and we accept the energy only if a uniform random number
    // is <= the acceptance fraction.
    do {
        // Generate an energy from the broken power law distribution
        energy = m_mc_brokenplaw.mc(emin, emax, time, ran);
        
        // Compute acceptance fraction
        acceptance_fraction = eval(energy) / m_mc_brokenplaw.eval(energy);
        
    } while (ran.uniform() > acceptance_fraction);
    
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
void GModelSpectralSmoothBrokenPlaw::read(const GXmlElement& xml)
{
    // Verify number of model parameters
    gammalib::xml_check_parnum(G_READ, xml, 6);

    // Get remaining XML parameters
    const GXmlElement* prefactor   = gammalib::xml_get_par(G_READ, xml, m_norm.name());
    const GXmlElement* index1      = gammalib::xml_get_par(G_READ, xml, m_index1.name());
    const GXmlElement* index2      = gammalib::xml_get_par(G_READ, xml, m_index2.name());
    const GXmlElement* pivot       = gammalib::xml_get_par(G_READ, xml, m_pivot.name());
    const GXmlElement* breakenergy = gammalib::xml_get_par(G_READ, xml, m_breakenergy.name());
    const GXmlElement* beta        = gammalib::xml_get_par(G_READ, xml, m_beta.name());
    
    // Read parameters
    m_norm.read(*prefactor);
    m_index1.read(*index1);
    m_index2.read(*index2);
    m_pivot.read(*pivot);
    m_breakenergy.read(*breakenergy);
    m_beta.read(*beta);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element.
 *
 * Writes the spectral information into an XML element.
 ***************************************************************************/
void GModelSpectralSmoothBrokenPlaw::write(GXmlElement& xml) const
{
    // Verify model type
    gammalib::xml_check_type(G_WRITE, xml, type());

    // Get XML parameters
    GXmlElement* norm        = gammalib::xml_need_par(G_WRITE, xml, m_norm.name());
    GXmlElement* index1      = gammalib::xml_need_par(G_WRITE, xml, m_index1.name());
    GXmlElement* index2      = gammalib::xml_need_par(G_WRITE, xml, m_index2.name());
    GXmlElement* pivot       = gammalib::xml_need_par(G_WRITE, xml, m_pivot.name());
    GXmlElement* breakenergy = gammalib::xml_need_par(G_WRITE, xml, m_breakenergy.name());
    GXmlElement* beta        = gammalib::xml_need_par(G_WRITE, xml, m_beta.name());
    
    // Write parameters
    m_norm.write(*norm);
    m_index1.write(*index1);
    m_pivot.write(*pivot);
    m_index2.write(*index2);
    m_breakenergy.write(*breakenergy);
    m_beta.write(*beta);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Print smooth broken powerlaw information
 *
 * @param[in] chatter Chattiness.
 * @return String containing model information.
 ***************************************************************************/
std::string GModelSpectralSmoothBrokenPlaw::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;
    
    // Continue only if chatter is not silent
    if (chatter != SILENT) {
        
        // Append header
        result.append("=== GModelSpectralSmoothBrokenPlaw ===");
        
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
void GModelSpectralSmoothBrokenPlaw::init_members(void)
{
    // Initialise model type
    m_type = "SmoothBrokenPowerLaw";

    // Initialise pre factor
    m_norm.clear();
    m_norm.name("Prefactor");
    m_norm.unit("ph/cm2/s/MeV");
    m_norm.scale(1.0);
    m_norm.value(1.0);          // default: 1.0
    m_norm.min(0.0);            // min:     0.0
    m_norm.free();
    m_norm.gradient(0.0);
    m_norm.has_grad(true);

    // Initialise index1
    m_index1.clear();
    m_index1.name("Index1");
    m_index1.scale(1.0);
    m_index1.value(-2.0);        // default: -2.0
    m_index1.range(-10.0,+10.0); // range:   [-10,+10]
    m_index1.free();
    m_index1.gradient(0.0);
    m_index1.has_grad(true);

    // Initialise index2
    m_index2.clear();
    m_index2.name("Index2");
    m_index2.scale(1.0);
    m_index2.value(-3.0);        // default: -2.0
    m_index2.range(-10.0,+10.0); // range:   [-10,+10]
    m_index2.free();
    m_index2.gradient(0.0);
    m_index2.has_grad(true);

    // Initialise pivot energy
    m_pivot.clear();
    m_pivot.name("PivotEnergy");
    m_pivot.unit("MeV");
    m_pivot.scale(1.0e5);
    m_pivot.value(1.0);          // default: 100 GeV
    m_pivot.fix();
    m_pivot.gradient(0.0);
    m_pivot.has_grad(true);

    // Initialise break energy
    m_breakenergy.clear();
    m_breakenergy.name("BreakEnergy");
    m_breakenergy.unit("MeV");
    m_breakenergy.scale(1.0e5);
    m_breakenergy.value(1.0);    // default: 100 GeV
    m_breakenergy.free();
    m_breakenergy.gradient(0.0);
    m_breakenergy.has_grad(true);

    // Initialize beta (break smoothness parameter)
    m_beta.clear();
    m_beta.name("BreakSmoothness");
    m_beta.scale(1.0);
    m_beta.value(1.0);
    m_beta.range(0.01,10.0);
    m_beta.free();
    m_beta.gradient(0.0);
    m_beta.has_grad(true);

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);
    m_pars.push_back(&m_index1);
    m_pars.push_back(&m_pivot);
    m_pars.push_back(&m_index2);
    m_pars.push_back(&m_breakenergy);
    m_pars.push_back(&m_beta);

    // Initialise eval cache
    m_last_energy.clear();
    m_last_index1          = 1.0e30;
    m_last_index2          = 1.0e30;
    m_last_pivot           = 1.0e30;
    m_last_breakenergy     = 1.0e30;
    m_last_beta            = 1.0e30;
    m_last_epivot_norm     = 1.0e30;
    m_last_ebreak_norm     = 1.0e30;
    m_last_log_epivot_norm = 1.0e30;
    m_last_log_ebreak_norm = 1.0e30;
    m_last_epivot_pow      = 1.0e30;
    m_last_ebreak_pow      = 1.0e30;

    // Initialise MC cache
    m_mc_prefactor   = 1.0e30;
    m_mc_index1      = 1.0e30;
    m_mc_index2      = 1.0e30;
    m_mc_pivot       = 1.0e30;
    m_mc_breakenergy = 1.0e30;
    m_mc_brokenplaw  = GModelSpectralBrokenPlaw();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Smooth broken power law.
 ***************************************************************************/
void GModelSpectralSmoothBrokenPlaw::copy_members(const GModelSpectralSmoothBrokenPlaw& model)
{
    // Copy members
    m_type        = model.m_type;
    m_norm        = model.m_norm;
    m_index1      = model.m_index1;
    m_index2      = model.m_index2;
    m_pivot       = model.m_pivot;
    m_breakenergy = model.m_breakenergy;
    m_beta        = model.m_beta;
    
    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);
    m_pars.push_back(&m_index1);
    m_pars.push_back(&m_pivot);
    m_pars.push_back(&m_index2);
    m_pars.push_back(&m_breakenergy);
    m_pars.push_back(&m_beta);
    
    // Copy eval cache
    m_last_energy          = model.m_last_energy;
    m_last_index1          = model.m_last_index1;
    m_last_index2          = model.m_last_index2;
    m_last_pivot           = model.m_last_pivot;
    m_last_breakenergy     = model.m_last_breakenergy;
    m_last_beta            = model.m_last_beta;
    m_last_epivot_norm     = model.m_last_epivot_norm;
    m_last_ebreak_norm     = model.m_last_ebreak_norm;
    m_last_log_epivot_norm = model.m_last_log_epivot_norm;
    m_last_log_ebreak_norm = model.m_last_log_ebreak_norm;
    m_last_epivot_pow      = model.m_last_epivot_pow;
    m_last_ebreak_pow      = model.m_last_ebreak_pow;
    
    // Copy MC cache
    m_mc_prefactor   = model.m_mc_prefactor;
    m_mc_index1      = model.m_mc_index1;
    m_mc_index2      = model.m_mc_index2;
    m_mc_pivot       = model.m_mc_pivot;
    m_mc_breakenergy = model.m_mc_breakenergy;
    m_mc_brokenplaw  = model.m_mc_brokenplaw;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralSmoothBrokenPlaw::free_members(void)
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
void GModelSpectralSmoothBrokenPlaw::update_eval_cache(const GEnergy& energy) const
{
    // Get parameter values
    double index1    = m_index1.value();
    double index2    = m_index2.value();
    double pivot_eng = pivot().MeV();
    double break_eng = breakenergy().MeV();
    double beta      = m_beta.value();
    
    // If the energy or one of the parameters index1, index2, breakenergy
    // energy, or beta has changed then recompute the cache
    if ((m_last_energy      != energy) ||
        (m_last_index1      != index1) ||
        (m_last_index2      != index2) ||
        (m_last_pivot       != pivot_eng) ||
        (m_last_breakenergy != break_eng) ||
        (m_last_beta        != beta)) {
        
        // Store actual energy and parameter values
        m_last_energy      = energy;
        m_last_index1      = index1;
        m_last_index2      = index2;
        m_last_pivot       = pivot_eng;
        m_last_breakenergy = break_eng;
        m_last_beta        = beta;
        
        // Compute and store value
        double eng             = energy.MeV();
        m_last_epivot_norm     = eng / m_last_pivot;
        m_last_ebreak_norm     = eng / m_last_breakenergy;
        m_last_log_epivot_norm = std::log(m_last_epivot_norm);
        m_last_log_ebreak_norm = std::log(m_last_ebreak_norm);
        
        m_last_epivot_pow = std::pow(m_last_epivot_norm,m_last_index1);
        m_last_ebreak_pow = std::pow(m_last_ebreak_norm,
                                     (m_last_index1-m_last_index2)/m_last_beta) ;
        
    } // endif: recomputation was required
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update Monte Carlo pre computation cache
 *
 * Updates the precomputation cache for Monte Carlo simulations.
 ***************************************************************************/
void GModelSpectralSmoothBrokenPlaw::update_mc_cache(void) const
{
    // Check if we need to update the cache
    if (prefactor()         != m_mc_prefactor ||
        index1()            != m_mc_index1    ||
        index2()            != m_mc_index2    ||
        pivot().MeV()       != m_mc_pivot     ||
        breakenergy().MeV() != m_mc_breakenergy) {

        // Set parameters for which Monte Carlo cache was computed
        m_mc_prefactor = prefactor();
        m_mc_index1      = index1();
        m_mc_index2      = index2();
        m_mc_pivot       = pivot().MeV();
        m_mc_breakenergy = breakenergy().MeV();

        // Compute prefactor at pivot energy
        double pre = prefactor() *
                     std::pow(breakenergy().MeV()/pivot().MeV(), index1());

        // Find out which index is harder. This is important since the smoothly
        // broken power law follows the hard index below the break energy and
        // the softer index above the break energy.
        double index1 = (m_mc_index1 > m_mc_index2) ? m_mc_index1 : m_mc_index2;
        double index2 = (m_mc_index1 > m_mc_index2) ? m_mc_index2 : m_mc_index1;

        // Set broken power law for Monte Carlo simulations
        m_mc_brokenplaw = GModelSpectralBrokenPlaw(pre,
                                                   index1,
                                                   breakenergy(),
                                                   index2);
        
    } // endif: Update was required
    
    // Return
    return;
}
