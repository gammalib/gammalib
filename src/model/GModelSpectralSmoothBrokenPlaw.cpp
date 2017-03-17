/***************************************************************************
 *               GModelSpectralSmoothBrokenPlaw.cpp -                      *
 *             Smoothly broken power law spectrum class                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Joshua Cardenzana                                *
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
 * @brief Broken power law spectrum class implementation
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
const GModelSpectralSmoothBrokenPlaw g_spectral_sblaw_seed2("SmoothBrokenPowerLaw",
                                                            "Prefactor",
                                                            "Index1",
                                                            "Scale",
                                                            "Index2",
                                                            "BreakValue",
                                                            "Beta");
const GModelSpectralRegistry  g_spectral_sblaw_registry1(&g_spectral_sblaw_seed1);
const GModelSpectralRegistry  g_spectral_sblaw_registry2(&g_spectral_sblaw_seed2);


/* __ Method name definitions ____________________________________________ */
#define G_MC   "GModelSpectralSmoothBrokenPlaw::mc(GEnergy&, GEnergy&, GTime&,"\
" GRan&)"
#define G_READ              "GModelSpectralSmoothBrokenPlaw::read(GXmlElement&)"
#define G_WRITE            "GModelSpectralSmoothBrokenPlaw::write(GXmlElement&)"

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
 * @param[in] type          Model type.
 * @param[in] prefactor     Name of prefactor parameter.
 * @param[in] index1        Name of index1 parameter.
 * @param[in] pivot         Name of pivot parameter.
 * @param[in] index2        Name of index2 parameter.
 * @param[in] breakenergy   Name of breakenergy parameter.
 * @param[in] beta          Name of beta parameter.
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
    m_breakenergy.name(breakenergy);
    m_index2.name(index2);
    m_beta.name(beta);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Parameter constructor
 *
 * @param[in] prefactor     Smoothly broken Plaw pre factor (ph/cm2/s/MeV).
 * @param[in] index1        Smoothly broken Plaw index1.
 * @param[in] pivot         Smoothly broken Plaw pivot energy
 * @param[in] index2        Smoothly broken Plaw index1.
 * @param[in] breakenergy   Break energy.
 * @param[in] beta          Break smoothness parameter
 *
 * Constructs a spectral power law using the model parameters
 *
 *     power law @p prefactor (ph/cm2/s/MeV)
 *     spectral @p index1
 *     Energy @p pivot energy
 *     spectral @p index2
 *     Energy @p breakenergy of spectral break
 *     spectral @p beta
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
 * Constructs a power law spectral model by extracting information from an
 * XML element. See the read() method for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelSpectralSmoothBrokenPlaw::GModelSpectralSmoothBrokenPlaw(
                                        const GXmlElement& xml) :
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
 *    S_{\rm E}(E | t) = k_0 \times \left \{
 *    \begin{eqnarray}
 *     \left( \frac{E}{E_b} \right)^{\gamma_1} & {\rm if\,\,} E < E_b \\
 *     \left( \frac{E}{E_b} \right)^{\gamma_2} & {\rm otherwise}
 *    \end{eqnarray}
 *    \right .
 * \f]
 *
 * where
 * \f$k_0\f$ is the normalization or prefactor,
 * \f$\gamma_1\f$ is the spectral index before the break,
 * \f$\gamma_2\f$ is the spectral index after the break, and
 * \f$E_b\f$ is the break energy.
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
 *    \frac{\delta S_{\rm E}(E | t)}{\delta \gamma_1} = \left \{
 *    \begin{eqnarray}
 *      S_{\rm E}(E | t) \, \ln(E/E_b) & {\rm if\,\,} E < E_b \\
 *      0 & {\rm otherwise}
 *    \end{eqnarray}
 *    \right .
 * \f]
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta \gamma_2} = \left \{
 *    \begin{eqnarray}
 *      0                              & {\rm if\,\,} E < E_b \\
 *      S_{\rm E}(E | t) \, \ln(E/E_b) & {\rm otherwise}
 *    \end{eqnarray}
 *    \right .
 * \f]
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta E_b} = k_0 \times \left \{
 *    \begin{eqnarray}
 *      -S_{\rm E}(E | t) \, \left( \frac{\gamma_1}{E_b} \right) & {\rm if\,\,} E < E_b \\
 *      -S_{\rm E}(E | t) \, \left( \frac{\gamma_2}{E_b} \right) & {\rm otherwise}
 *    \end{eqnarray}
 *    \right .
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
                ? value * m_index1.scale() / (1.0+m_last_ebreak_pow) *
                  (m_last_log_epivot_norm+(m_last_log_epivot_norm-m_last_log_ebreak_norm)*
                  m_last_ebreak_pow)
                : 0.0;
        double g_index2 = (m_index2.is_free())
                ? value * m_index2.scale() / (1.0+m_last_ebreak_pow) *
                  m_last_log_ebreak_norm * m_last_ebreak_pow
                : 0.0;
        // Compute pivot and break energy value gradients
        double g_pivot = (m_pivot.is_free())
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
 * @exception GException::erange_invalid
 *            Energy range is invalid (emin < emax required).
 *
 * Returns Monte Carlo energy by randomly drawing from a smoothly broken power law.
 ***************************************************************************/
GEnergy GModelSpectralSmoothBrokenPlaw::mc(const GEnergy& emin,
                                           const GEnergy& emax,
                                           const GTime&   time,
                                           GRan&          ran) const
{
    // Throw exception if energy range is not valid
    if (emin >= emax) {
        throw GException::erange_invalid(G_MC, emin.MeV(), emax.MeV(),
                            "Minimum energy < maximum energy required.");
    }
    
    // Allocate energy
    GEnergy energy;
    
    // Update cache
    update_mc_cache(emin, emax);
    
    // Initialise energy
    double eng(0.0);
    
    // Initialse acceptance fraction
    double acceptance_fraction(0.0);
    
    // Use rejection method to draw a random energy. We first draw
    // analytically from a power law, and then compare the power law
    // at the drawn energy to the curved function. This
    // gives an acceptance fraction, and we accept the energy only if
    // a uniform random number is <= the acceptance fraction.
    do {
        
        double e_norm(0.0), plaw(0.0);
        
        // Get uniform random number
        double u = ran.uniform();
        
        // Case 1: indices are the same (things are simple)
        if (m_mc_exponentH == m_mc_exponentS) {
            // Case 1A: Both indices are -1
            if (m_mc_exponentH == 0.0) {
                eng = std::exp(u * m_mc_norm + std::log(m_mc_emin));
                
            }
            // Case 1B: Indices are not equal to -1
            else {
                eng = std::exp(std::log(u * m_mc_norm +
                      std::pow(m_mc_emin,m_mc_exponentH)) / m_mc_exponentH);
            }
            // There is no differentiation between the powerlaws when the
            // exponents are the same.
            e_norm = eng / m_breakenergy.value();
            plaw = m_mc_plaw_prefactor * std::pow(e_norm, m_mc_exponentH-1.0);
            
        }
        // Case 2: indices are different (things are more complicated)
        else {
        
            // Case 2A: Random number suggests energy less than breakenergy
            if (u <= (m_mc_pow_ewidth_low/m_mc_norm)) {
                
                // Case 2A-1: Hard index is not -1
                if (m_mc_exponentH != 0.0) {
                    eng = std::exp(std::log(u * m_mc_norm +
                          std::pow(m_mc_emin, m_mc_exponentH)) / m_mc_exponentH);
                }
                
                // Case 2A-2: Hard index is -1
                else {
                    eng = std::exp(u * m_mc_norm + std::log(m_mc_emin));
                }
                
                // Compute powerlaw at given energy
                e_norm = eng / m_breakenergy.value();
                plaw   = m_mc_plaw_prefactor * std::pow(e_norm, m_mc_exponentH-1.0);
                
            }
            // Case 2B: Random number suggests energy above breakenergy
            else {
                // Case 2B-1: Soft index is not -1
                if (m_mc_exponentS != 0.0) {
                    eng = std::exp(std::log(u * m_mc_norm +
                          std::pow(breakenergy().MeV(), m_mc_exponentS) -
                          m_mc_pow_ewidth_low) / m_mc_exponentS) ;
                }
                // Case 2B-2: Soft index is -1
                else {
                    eng = std::exp(u * m_mc_norm +
                                   std::log(m_breakenergy.value()) -
                                   m_mc_pow_ewidth_low);
                }
                
                // Compute powerlaw at given energy
                e_norm = eng / m_breakenergy.value();
                plaw   = m_mc_plaw_prefactor * std::pow(e_norm, m_mc_exponentS-1.0);
                
            }
        }
        
        // Compute logparabola at given energy
        double sb_plaw = prefactor() * std::pow(eng/m_pivot.value(),index1()) *
            std::pow(1.0 + std::pow(e_norm,(index1()-index2())/beta()), -beta());
        
        // Compute acceptance fraction
        acceptance_fraction = sb_plaw / plaw;
        
    } while (ran.uniform() > acceptance_fraction);
    
    // Set energy
    energy.MeV(eng);
    
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
    // Get remaining XML parameters
    const GXmlElement* prefactor   = gammalib::xml_get_par(G_READ, xml, m_norm.name());
    const GXmlElement* index1      = gammalib::xml_get_par(G_READ, xml, m_index1.name());
    const GXmlElement* pivot       = gammalib::xml_get_par(G_READ, xml, m_pivot.name());
    const GXmlElement* index2      = gammalib::xml_get_par(G_READ, xml, m_index2.name());
    const GXmlElement* breakenergy = gammalib::xml_get_par(G_READ, xml, m_breakenergy.name());
    const GXmlElement* beta        = gammalib::xml_get_par(G_READ, xml, m_beta.name());
    
    // Read parameters
    m_norm.read(*prefactor);
    m_index1.read(*index1);
    m_pivot.read(*pivot);
    m_index2.read(*index2);
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
 * @exception GException::model_invalid_spectral
 *            Existing XML element is not of the expected type.
 *
 * Writes the spectral information into an XML element.
 ***************************************************************************/
void GModelSpectralSmoothBrokenPlaw::write(GXmlElement& xml) const
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
 * @param[in] chatter Chattiness (defaults to NORMAL).
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
    m_pivot.clear();
    m_pivot.name("PivotEnergy");
    m_pivot.unit("MeV");
    m_pivot.scale(1.0);
    m_pivot.value(1.0e5);  // default: 100 GeV
    m_pivot.fix();
    m_pivot.gradient(0.0);
    m_pivot.has_grad(true);
    
    // Initialise break energy
    m_breakenergy.clear();
    m_breakenergy.name("BreakEnergy");
    m_breakenergy.unit("MeV");
    m_breakenergy.scale(1.0);
    m_breakenergy.value(100.0);  // default: 100
    m_breakenergy.fix();
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
    m_last_index1      = 1.0e30;
    m_last_index2      = 1.0e30;
    m_last_breakenergy = 1.0e30;
    m_last_beta        = 1.0e30;
    m_last_epivot_norm = 1.0e30;
    m_last_ebreak_norm = 1.0e30;
    m_last_log_epivot_norm = 1.0e30;
    m_last_log_ebreak_norm = 1.0e30;
    m_last_epivot_pow  = 1.0e30;
    m_last_ebreak_pow  = 1.0e30;
    
    // Initialise MC cache
    m_mc_emin           = 0.0;
    m_mc_emax           = 0.0;
    m_mc_plaw_prefactor = 0.0;
    m_mc_exponentS      = 0.0;
    m_mc_exponentH      = 0.0;
    m_mc_pow_ewidth_low = 0.0;
    m_mc_norm           = 0.0;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model GModelSpectralSmoothBrokenPlaw members which should be copied.
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
    m_pars.push_back(&m_index2);
    m_pars.push_back(&m_pivot);
    m_pars.push_back(&m_breakenergy);
    m_pars.push_back(&m_beta);
    
    // Copy eval cache
    m_last_energy      = model.m_last_energy;
    m_last_index1      = model.m_last_index1;
    m_last_index2      = model.m_last_index2;
    m_last_pivot       = model.m_last_pivot;
    m_last_breakenergy = model.m_last_breakenergy;
    m_last_beta        = model.m_last_beta;
    m_last_epivot_norm = model.m_last_epivot_norm;
    m_last_ebreak_norm = model.m_last_ebreak_norm;
    m_last_log_epivot_norm = model.m_last_log_epivot_norm;
    m_last_log_ebreak_norm = model.m_last_log_ebreak_norm;
    m_last_epivot_pow  = model.m_last_epivot_pow;
    m_last_ebreak_pow  = model.m_last_ebreak_pow;
    
    // Copy MC cache
    m_mc_emin           = model.m_mc_emin;
    m_mc_emax           = model.m_mc_emax;
    m_mc_plaw_prefactor = model.m_mc_plaw_prefactor;
    m_mc_exponentS      = model.m_mc_exponentS;
    m_mc_exponentH      = model.m_mc_exponentH;
    m_mc_pow_ewidth_low = model.m_mc_pow_ewidth_low;
    m_mc_norm           = model.m_mc_norm;
    
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
    // Get parameter values (takes 2 multiplications which are difficult
    // to avoid)
    double index1      = m_index1.value();
    double index2      = m_index2.value();
    double pivot       = m_pivot.value();
    double breakenergy = m_breakenergy.value();
    double beta        = m_beta.value();
    
    // If the energy or one of the parameters index1, index2, breakenergy
    // energy, or beta has changed then recompute the cache
    if ((m_last_energy       != energy) ||
        (m_last_index1       != index1) ||
        (m_last_index2       != index2) ||
        (m_last_pivot        != pivot)  ||
        (m_last_breakenergy  != breakenergy) ||
        (m_last_beta         != beta)) {
        
        // Store actual energy and parameter values
        m_last_energy       = energy;
        m_last_index1       = index1;
        m_last_index2       = index2;
        m_last_pivot        = pivot;
        m_last_breakenergy  = breakenergy;
        m_last_beta         = beta;
        
        // Compute and store value
        double eng        = energy.MeV();
        m_last_epivot_norm = eng / m_last_pivot;
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
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 *
 * Updates the precomputation cache for Monte Carlo simulations.
 ***************************************************************************/
void GModelSpectralSmoothBrokenPlaw::update_mc_cache(const GEnergy& emin,
                                                     const GEnergy& emax) const
{
    // Update the eval cache
    update_eval_cache(emin);
    
    // Check if we need to update the cache
    if (emin.MeV() != m_mc_emin || emax.MeV() != m_mc_emax) {
        
        // Store new energy interval
        m_mc_emin = emin.MeV();
        m_mc_emax = emax.MeV();
        
        /* Predefine some variables */
        
        // Prefactor for comparison power law functions
        m_mc_plaw_prefactor = prefactor() * std::pow(m_last_breakenergy/m_last_pivot,
                                                     m_last_index1);
        
        // Find out which index is harder. This is important since the
        // smoothly broken power law follows the hard index below the break
        // energy and the softer index above the break energy.
        
        // CASE: index1 is harder
        if (index1() > index2()) {
            m_mc_exponentH = m_last_index1 + 1.0 ;   // Hard index + 1
            m_mc_exponentS = m_last_index2 + 1.0 ;   // Soft index + 1
            
        }
        // CASE: index2 is harder, or index1==index2
        else {
            m_mc_exponentH = m_last_index2 + 1.0 ;   // Hard index + 1
            m_mc_exponentS = m_last_index1 + 1.0 ;   // Soft index + 1
        }
        
        // Now handle the cases where the hard exponent is -1
        if (m_mc_exponentH == 0.0) {
            
            // CASE: both exponents are -1
            if (m_mc_exponentS == 0.0) {
                m_mc_pow_ewidth_low = 0.0 ;
                m_mc_norm = std::log(m_mc_emax) - std::log(m_mc_emin);
            }
            // CASE: Only harder exponent is -1
            else {
                m_mc_pow_ewidth_low = std::log(m_last_breakenergy/m_mc_emin);
                m_mc_norm = m_mc_pow_ewidth_low + std::pow(m_mc_emax,m_mc_exponentS) -
                            std::pow(m_last_breakenergy,m_mc_exponentS);
            }
        }
        
        // CASE: Only softer exponent is -1
        else if (m_mc_exponentS == 0.0) {
            m_mc_pow_ewidth_low = std::pow(m_last_breakenergy,m_mc_exponentH) -
                                  std::pow(m_mc_emin, m_mc_exponentH) ;
            m_mc_norm = m_mc_pow_ewidth_low + std::log(m_mc_emax) -
                        std::log(m_last_breakenergy) ;
        }
        
        // CASE: Neither exponent is -1.
        else {
            // Get the width of the energy range below the break energy
            m_mc_pow_ewidth_low = std::pow(m_last_breakenergy, m_mc_exponentH) -
                                  std::pow(m_mc_emin, m_mc_exponentH);
            
            // Get the normalization term
            m_mc_norm = m_mc_pow_ewidth_low + (std::pow(m_mc_emax,m_mc_exponentS) -
                                               std::pow(m_last_breakenergy,m_mc_exponentS));
        }
    } // endif: Update was required
    
    // Return
    return;
}
