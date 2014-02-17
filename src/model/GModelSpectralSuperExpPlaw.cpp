/***************************************************************************
 *   GModelSpectralSuperExpPlaw.cpp - Super exp. cut off power law model   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Michael Mayer                                    *
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
 * @file GModelSpectralSuperExpPlaw.cpp
 * @brief Super exponential cut off power law spectral class implementation
 * @author Michael Mayer
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GIntegral.hpp"
#include "GModelSpectralSuperExpPlaw.hpp"
#include "GModelSpectralRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpectralSuperExpPlaw  g_spectral_supexpplaw_seed;
const GModelSpectralRegistry      g_spectral_supexpplaw_registry(&g_spectral_supexpplaw_seed);

/* __ Method name definitions ____________________________________________ */
#define G_MC     "GModelSpectralSuperExpPlaw::mc(GEnergy&, GEnergy&, GTime&,"\
                                                                    " GRan&)"
#define G_READ               "GModelSpectralSuperExpPlaw::read(GXmlElement&)"
#define G_WRITE             "GModelSpectralSuperExpPlaw::write(GXmlElement&)"

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
GModelSpectralSuperExpPlaw::GModelSpectralSuperExpPlaw(void) : GModelSpectral()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] prefactor Pre factor normalization (ph/cm2/s/MeV).
 * @param[in] index1 Power law index.
 * @param[in] pivot Pivot energy.
 * @param[in] cutoff Cut off energy.
 * @param[in] index2 index of cut off
 *
 * Construct an exponentially cut off power law from
 * - a prefactor value (in units of ph/cm2/s/MeV)
 * - a spectral index,
 * - a pivot energy,
 * - a cut off energy, and
 * - a cut off exponent
 ***************************************************************************/
GModelSpectralSuperExpPlaw::GModelSpectralSuperExpPlaw(const double&  prefactor,
                                                       const double&  index1,
                                                       const GEnergy& pivot,
                                                       const GEnergy& cutoff,
                                                       const double&  index2) :
                            GModelSpectral()
{
    // Initialise members
    init_members();

    // Set parameters
    m_norm.value(prefactor);
    m_index1.value(index1);
    m_pivot.value(pivot.MeV()); // Internally stored in MeV
    m_ecut.value(cutoff.MeV()); // Internally stored in MeV
    m_index2.value(index2);

    // Autoscale parameters
    autoscale();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs an exponentially cut off power law spectral model by extracting
 * information from an XML element. See the read() method for more
 * information about the expected structure of the XML element.
 ***************************************************************************/
GModelSpectralSuperExpPlaw::GModelSpectralSuperExpPlaw(const GXmlElement& xml) :
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
 * @param[in] model Super exponentially cut off power law model.
 ***************************************************************************/
GModelSpectralSuperExpPlaw::GModelSpectralSuperExpPlaw(const GModelSpectralSuperExpPlaw& model) :
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
GModelSpectralSuperExpPlaw::~GModelSpectralSuperExpPlaw(void)
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
 * @param[in] model Super exponentially cut off power law model.
 * @return Super exponentially cut off power law model.
 ***************************************************************************/
GModelSpectralSuperExpPlaw& GModelSpectralSuperExpPlaw::operator=(const GModelSpectralSuperExpPlaw& model)
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
 * @brief Clear exponentially cut off power law model
 ***************************************************************************/
void GModelSpectralSuperExpPlaw::clear(void)
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
 * @brief Clone super exponentially cut off power law model
 *
 * @return Pointer to deep copy of super exponentially cut off power law model.
 ***************************************************************************/
GModelSpectralSuperExpPlaw* GModelSpectralSuperExpPlaw::clone(void) const
{
    // Clone exponentially cut off power law model
    return new GModelSpectralSuperExpPlaw(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @return Model value.
 *
 * Evaluates
 *
 * \f[
 *    S_{\rm E}(E | t) = {\tt m\_norm}
 *    \left( \frac{E}{\tt m\_pivot} \right)^{\tt m\_index1}
 *    \exp \left( \frac{-E}{\tt m\_ecut}^{\tt m_index2} \right)
 * \f]
 *
 * where
 * - \f${\tt m\_norm}\f$ is the normalization or prefactor,
 * - \f${\tt m\_pivot}\f$ is the pivot energy,
 * - \f${\tt m\_index1}\f$ is the spectral index
 * - \f${\tt m\_ecut}\f$ is the cut off energy, and
 * - \f${\tt m\_index2}\f$ is the index of the cut off
 *
 * @todo The method expects that pivot!=0 and ecut!=0. Otherwise Inf or NaN
 *       may result. We should add a test that prevents using invalid
 *       values.
 ***************************************************************************/
double GModelSpectralSuperExpPlaw::eval(const GEnergy& srcEng,
                                        const GTime&   srcTime) const
{
    // Update the evaluation cache
    update_eval_cache(srcEng);

    // Compute function value
    double value = m_norm.value() * m_last_power;

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GModelSpectralSuperExpPlaw::eval";
        std::cout << "(srcEng=" << srcEng;
        std::cout << ", srcTime=" << srcTime << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", power=" << m_last_power;
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
 * @return Model value.
 *
 * Evaluates
 *
 * \f[
 *    S_{\rm E}(E | t) = {\tt m\_norm}
 *    \left( \frac{E}{\tt m\_pivot} \right)^{\tt m\_index1}
 *    \exp \left( \frac{-E}{\tt m\_ecut}^{\tt m\_index2} \right)
 * \f]
 *
 * where
 * - \f${\tt m\_norm}\f$ is the normalization or prefactor,
 * - \f${\tt m\_index1}\f$ is the spectral index,
 * - \f${\tt m\_ecut}\f$ is the cut off energy
 * - \f${\tt m\_pivot}\f$ is the pivot energy, and
 * - \f${\tt m\_index2}\f$ is the index of the cut off
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
 *    \frac{\delta S_{\rm E}(E | t)}{\delta {\tt m\_index1}} =
 *      S_{\rm E}(E | t) \, \ln(E/{\tt m_pivot})
 * \f]
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta {\tt m\_ecut}} =
 *      S_{\rm E}(E | t) \, \left( \frac{E}{{\tt m\_ecut}^{\tt m\_index2}} \right) \, \left({\frac{\tt m\_index2}}{\tt m\_ecut}\right)
 * \f]
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta {\tt m\_index2}} =
 *      - S_{\rm E}(E | t) \, \ln \left(\frac{E}{{\tt m\_ecut}\right)\, \left(\frac{E}{{\tt m\_ecut}^{\tt m\_index2}} \right)
 * \f]
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta {\tt m\_pivot}} =
 *      -S_{\rm E}(E | t) \,
 *      \left( \frac{{\tt m\_index}}{{\tt m\_pivot}} \right)
 * \f]
 *
 * @todo The method expects that pivot!=0 and ecut!=0. Otherwise Inf or NaN
 *       may result. We should add a test that prevents using invalid
 *       values.
 ***************************************************************************/
double GModelSpectralSuperExpPlaw::eval_gradients(const GEnergy& srcEng,
                                                  const GTime&   srcTime)
{
    // Update the evaluation cache
    update_eval_cache(srcEng);

    // Compute function value
    double value = m_norm.value() * m_last_power;

    // Compute partial derivatives with respect to the parameter factor
    // values. The partial derivatives with respect to the parameter
    // values are obtained by division by the scale factor.
    double g_norm  = (m_norm.is_free())
                     ? m_norm.scale() * m_last_power : 0.0;
    double g_index1 = (m_index1.is_free())
                     ? value * m_index1.scale() * std::log(m_last_e_norm) : 0.0;
    double g_ecut  = (m_ecut.is_free())
                     ? value * m_last_index2 * m_last_exponent / m_ecut.factor_value() : 0.0;
    double g_pivot = (m_pivot.is_free())
                     ? -value * m_last_index1 / m_pivot.factor_value() : 0.0;
    double g_index2 = (m_index2.is_free())
                      ? - value * m_index2.scale() * std::log(m_last_e_cut) * m_last_exponent : 0.0;

    // Set gradients
    m_norm.factor_gradient(g_norm);
    m_index1.factor_gradient(g_index1);
    m_ecut.factor_gradient(g_ecut);
    m_pivot.factor_gradient(g_pivot);
    m_index2.factor_gradient(g_index2);

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GModelSpectralSuperExpPlaw::eval_gradients";
        std::cout << "(srcEng=" << srcEng;
        std::cout << ", srcTime=" << srcTime << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", e_norm=" << m_last_e_norm;
        std::cout << ", e_cut=" << m_last_e_cut;
        std::cout << ", power=" << m_last_power;
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
double GModelSpectralSuperExpPlaw::flux(const GEnergy& emin,
                                        const GEnergy& emax) const
{
    // Initialise flux
    double flux = 0.0;
    
    // Compute only if integration range is valid
    if (emin < emax) {

        // Setup integration kernel
        flux_kernel integrand(m_norm.value(),  m_index1.value(),
                              m_pivot.value(), m_ecut.value(),m_index2.value());
        GIntegral integral(&integrand);

        // Get integration boundaries in MeV
        double e_min = emin.MeV();
        double e_max = emax.MeV();

        // Perform integration
        flux = integral.romb(e_min, e_max);

    } // endif: integration range was valid

    // Return
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
double GModelSpectralSuperExpPlaw::eflux(const GEnergy& emin,
                                         const GEnergy& emax) const
{
    // Initialise flux
    double eflux = 0.0;
    
    // Compute only if integration range is valid
    if (emin < emax) {

        // Setup integration kernel
        eflux_kernel integrand(m_norm.value(),  m_index1.value(),
                               m_pivot.value(), m_ecut.value(), m_index2.value());
        GIntegral integral(&integrand);

        // Get integration boundaries in MeV
        double e_min = emin.MeV();
        double e_max = emax.MeV();

        // Perform integration
        eflux = integral.romb(e_min, e_max);

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
 * Simulates a random energy in the interval [emin, emax] for an
 * exponentially cut off power law. The simulation is done using a rejection
 * method. First, a random energy within [emin, emax] is drawn from an
 * power law distribution. Then the energy is accepted or rejected based
 * on an acceptance fraction that is computed from the exponential cut off.
 ***************************************************************************/
GEnergy GModelSpectralSuperExpPlaw::mc(const GEnergy& emin,
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

    // Initialise energy
    double eng;
    
    // Initialse acceptance fraction
    double acceptance_fraction;
    double inv_ecut = 1.0 / m_ecut.value();

    // Use rejection method to draw a random energy. We first draw
    // analytically from a power law, and then compare the power law
    // at the drawn energy to the exponentially cut off function. This
    // gives an acceptance fraction, and we accept the energy only if
    // a uniform random number is <= the acceptance fraction.
    do {

        // Get uniform random number
        double u = ran.uniform();

        // Case A: Index is not -1
        if (index1() != -1.0) {
            if (u > 0.0) {
                eng = std::exp(std::log(u * m_mc_pow_ewidth + m_mc_pow_emin) /
                               m_mc_exponent);
            }
            else {
                eng = 0.0;
            }
        }

        // Case B: Index is -1
        else {
            eng = std::exp(u * m_mc_pow_ewidth + m_mc_pow_emin);
        }
        
        // Compute acceptance fraction
        acceptance_fraction = std::exp(- std::pow(eng * inv_ecut, m_index2.value()));

    } while (ran.uniform() > acceptance_fraction);
    
    // Set energy
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
 *     <spectrum type="PLSuperExpCutoff">
 *       <parameter name="Prefactor" scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Index1"    scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Cutoff"    scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Scale"     scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Index2"    scale=".." value=".." min=".." max=".." free=".."/>
 *     </spectrum>
 ***************************************************************************/
void GModelSpectralSuperExpPlaw::read(const GXmlElement& xml)
{
    // Verify that XML element has exactly 5 parameters
    if (xml.elements() != 5 || xml.elements("parameter") != 5) {
        throw GException::model_invalid_parnum(G_READ, xml,
              "Super Exponential Power law model requires exactly 5 parameters.");
    }

    // Extract model parameters
    int npar[] = {0, 0, 0, 0, 0};
    for (int i = 0; i < 5; ++i) {

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

        // Handle cutoff
        else if (par->attribute("name") == "Cutoff") {
            m_ecut.read(*par);
            npar[2]++;
        }

        // Handle pivot energy
        else if (par->attribute("name") == "Scale") {
            m_pivot.read(*par);
            npar[3]++;
        }

        // Handle index2
        else if (par->attribute("name") == "Index2") {
            m_index2.read(*par);
            npar[4]++;
        }

    } // endfor: looped over all parameters

    // Verify that all parameters were found
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1 || npar[3] != 1 || npar[4] != 1) {
        throw GException::model_invalid_parnames(G_READ, xml,
              "Require \"Prefactor\", \"Index1\", \"Cutoff\", \"Scale\" and \"Index2\""
              " parameters.");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * @exception GException::model_invalid_spectral
 *            Existing XML element is not of type "ExpCutoff"
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters or nodes found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Writes the spectral information into an XML element. The format of the XML
 * element is
 *
 *     <spectrum type="PLSuperExpCutoff">
 *       <parameter name="Prefactor" scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Index1"     scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Cutoff"    scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Scale"     scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Index2"     scale=".." value=".." min=".." max=".." free=".."/>
 *     </spectrum>
 ***************************************************************************/
void GModelSpectralSuperExpPlaw::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "") {
        xml.attribute("type", "PLSuperExpCutoff");
    }

    // Verify model type
    if (xml.attribute("type") != "PLSuperExpCutoff") {
        throw GException::model_invalid_spectral(G_WRITE, xml.attribute("type"),
              "Spectral model is not of type \"PLSuperExpCutoff\".");
    }

    // If XML element has 0 nodes then append 5 parameter nodes
    if (xml.elements() == 0) {
        xml.append(GXmlElement("parameter name=\"Prefactor\""));
        xml.append(GXmlElement("parameter name=\"Index1\""));
        xml.append(GXmlElement("parameter name=\"Cutoff\""));
        xml.append(GXmlElement("parameter name=\"Scale\""));
        xml.append(GXmlElement("parameter name=\"Index2\""));
    }

    // Verify that XML element has exactly 5 parameters
    if (xml.elements() != 5 || xml.elements("parameter") != 5) {
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Super Exponetial Power law model requires exactly 5 parameters.");
    }

    // Set or update model parameter attributes
    int npar[] = {0, 0, 0, 0, 0};
    for (int i = 0; i < 5; ++i) {

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

        // Handle index
        else if (par->attribute("name") == "Cutoff") {
            npar[2]++;
            m_ecut.write(*par);
        }

        // Handle pivot energy
        else if (par->attribute("name") == "Scale") {
            m_pivot.write(*par);
            npar[3]++;
        }

        // Handle index2
        else if (par->attribute("name") == "Index2") {
            npar[4]++;
            m_index2.write(*par);
        }

    } // endfor: looped over all parameters

    // Check of all required parameters are present
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1 || npar[3] != 1 || npar[4] != 1) {
        throw GException::model_invalid_parnames(G_WRITE, xml,
              "Require \"Prefactor\", \"Index1\", \"Cutoff\", \"Scale\" and \"Index2\""
              " parameters.");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print model information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing model information.
 ***************************************************************************/
std::string GModelSpectralSuperExpPlaw::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpectralSuperExpPlaw ===");

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
void GModelSpectralSuperExpPlaw::init_members(void)
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
    m_norm.has_grad(true);

    // Initialise powerlaw index
    m_index1.clear();
    m_index1.name("Index1");
    m_index1.scale(1.0);
    m_index1.value(-2.0);        // default: -2.0
    m_index1.range(-10.0,+10.0); // range:   [-10,+10]
    m_index1.free();
    m_index1.gradient(0.0);
    m_index1.has_grad(true);

    // Initialise cut off energy
    m_ecut.clear();
    m_ecut.name("Cutoff");
    m_ecut.unit("MeV");
    m_ecut.scale(1.0);
    m_ecut.value(1000.0);       // default: 1000.0
    m_ecut.min(0.1);            // min:     0.1
    m_ecut.free();
    m_ecut.gradient(0.0);
    m_ecut.has_grad(true);

    // Initialise pivot energy
    m_pivot.clear();
    m_pivot.name("PivotEnergy");
    m_pivot.unit("MeV");
    m_pivot.scale(1.0);
    m_pivot.value(100.0);
    m_pivot.fix();
    m_pivot.gradient(0.0);
    m_pivot.has_grad(true);

    // Initialise index2
    m_index2.clear();
	m_index2.name("Index2");
	m_index2.scale(1.0);
	m_index2.value(1.0);        // default: -2.0
	m_index2.range(0.0,5.0); // range:   [-10,+10]
	m_index2.free();
	m_index2.gradient(0.0);
	m_index2.has_grad(true);

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);
    m_pars.push_back(&m_index1);
    m_pars.push_back(&m_ecut);
    m_pars.push_back(&m_pivot);
    m_pars.push_back(&m_index2);

    // Initialise eval cache
    m_last_energy.clear();
    m_last_index1  = 1.0e30;
    m_last_ecut   = 1.0e30;
    m_last_pivot  = 1.0e30;
    m_last_index2 = 1.0e30;
    m_last_e_norm = 0.0;
    m_last_e_cut  = 0.0;
    m_last_power  = 0.0;
    m_last_exponent = 0.0;

    // Initialise MC cache
    m_mc_emin       = 0.0;
    m_mc_emax       = 0.0;
    m_mc_exponent   = 0.0;
    m_mc_pow_emin   = 0.0;
    m_mc_pow_ewidth = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Exponential cut off power law model.
 ***************************************************************************/
void GModelSpectralSuperExpPlaw::copy_members(const GModelSpectralSuperExpPlaw& model)
{
    // Copy members
    m_norm  = model.m_norm;
    m_index1 = model.m_index1;
    m_ecut  = model.m_ecut;
    m_pivot = model.m_pivot;
    m_index2 = model.m_index2;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);
    m_pars.push_back(&m_index1);
    m_pars.push_back(&m_ecut);
    m_pars.push_back(&m_pivot);
    m_pars.push_back(&m_index2);

    // Copy eval cache
    m_last_energy = model.m_last_energy;
    m_last_index1  = model.m_last_index1;
    m_last_ecut   = model.m_last_ecut;
    m_last_pivot  = model.m_last_pivot;
    m_last_index1  = model.m_last_index1;
    m_last_e_norm = model.m_last_e_norm;
    m_last_e_cut  = model.m_last_e_cut;
    m_last_power  = model.m_last_power;
    m_last_exponent = model.m_last_exponent;

    // Copy MC cache
    m_mc_emin       = model.m_mc_emin;
    m_mc_emax       = model.m_mc_emax;
    m_mc_exponent   = model.m_mc_exponent;
    m_mc_pow_emin   = model.m_mc_pow_emin;
    m_mc_pow_ewidth = model.m_mc_pow_ewidth;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralSuperExpPlaw::free_members(void)
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
void GModelSpectralSuperExpPlaw::update_eval_cache(const GEnergy& energy) const
{
    // Get parameter values (takes 3 multiplications which are difficult
    // to avoid)
    double index1 = m_index1.value();
    double ecut  = m_ecut.value();
    double pivot = m_pivot.value();
    double index2 = m_index2.value();
    
    // If the energy or one of the parameters index1, index2, cut-off or pivot
    // energy has changed then recompute the cache
    if ((m_last_energy != energy) ||
        (m_last_index1 != index1) ||
        (m_last_ecut   != ecut)   ||
        (m_last_pivot  != pivot)  ||
        (m_last_index2 != index2)) {

        // Store actual energy and parameter values
        m_last_energy = energy;
        m_last_index1 = index1;
        m_last_ecut   = ecut;
        m_last_pivot  = pivot;
        m_last_index2 = index2;

        // Compute and store value
        double eng      = energy.MeV();
        m_last_e_norm   = eng / m_last_pivot;
        m_last_e_cut    = eng / m_last_ecut;
        m_last_exponent = std::pow(m_last_e_cut, m_last_index2);
        m_last_power    = std::pow(m_last_e_norm, m_last_index1) *
                          std::exp(-m_last_exponent);

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
void GModelSpectralSuperExpPlaw::update_mc_cache(const GEnergy& emin,
                                            const GEnergy& emax) const

{
    // Case A: Index is not -1
    if (index1() != -1.0) {

        // Change in energy boundaries?
        if (emin.MeV() != m_mc_emin || emax.MeV() != m_mc_emax) {
            m_mc_emin       = emin.MeV();
            m_mc_emax       = emax.MeV();
            m_mc_exponent   = index1() + 1.0;
            m_mc_pow_emin   = std::pow(m_mc_emin, m_mc_exponent);
            m_mc_pow_ewidth = std::pow(m_mc_emax, m_mc_exponent) - m_mc_pow_emin;
        }

    }

    // Case B: Index is -1
    else {

        // Change in energy boundaries?
        if (emin.MeV() != m_mc_emin || emax.MeV() != m_mc_emax) {
            m_mc_emin       = emin.MeV();
            m_mc_emax       = emax.MeV();
            m_mc_exponent   = 0.0;
            m_mc_pow_emin   = std::log(m_mc_emin);
            m_mc_pow_ewidth = std::log(m_mc_emax) - m_mc_pow_emin;
        }

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Kernel for photon flux integration
 *
 * @param[in] energy Energy (MeV).
 ***************************************************************************/
double GModelSpectralSuperExpPlaw::flux_kernel::eval(const double& energy)
{
    // Evaluate function value
    double e_norm   = energy * m_inv_pivot;
    double e_cut    = energy * m_inv_ecut;
    double exponent = std::pow(e_cut, m_index2);
    double power    = std::pow(e_norm, m_index1) * std::exp(-exponent);
    double value    = m_norm * power;

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Kernel for energy flux integration
 *
 * @param[in] energy Energy (MeV).
 ***************************************************************************/
double GModelSpectralSuperExpPlaw::eflux_kernel::eval(const double& energy)
{
    // Evaluate function value
    double e_norm   = energy * m_inv_pivot;
    double e_cut    = energy * m_inv_ecut;
    double exponent = std::pow(e_cut, m_index2);
    double power    = std::pow(e_norm, m_index1) * std::exp(-exponent);
    double value    = m_norm * power * energy;

    // Return value
    return value;
}
