/***************************************************************************
 *     GModelSpectralExpPlaw.cpp - Exponential cut off power law model     *
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
 * @file GModelSpectralExpPlaw2.cpp
 * @brief Exponential cut off power law spectral class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GRan.hpp"
#include "GIntegral.hpp"
#include "GModelSpectralExpPlaw2.hpp"
#include "GModelSpectralRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpectralExpPlaw2 g_spectral_expplaw2_seed;
const GModelSpectralRegistry g_spectral_expplaw2_registry(&g_spectral_expplaw2_seed);

/* __ Method name definitions ____________________________________________ */
#define G_MC  "GModelSpectralExpPlaw2::mc(GEnergy&, GEnergy&, GTime&, GRan&)"
#define G_READ                   "GModelSpectralExpPlaw2::read(GXmlElement&)"
#define G_WRITE                 "GModelSpectralExpPlaw2::write(GXmlElement&)"

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
GModelSpectralExpPlaw2::GModelSpectralExpPlaw2(void) : GModelSpectral()
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
 * @param[in] index Power law index.
 * @param[in] pivot Pivot energy.
 * @param[in] lambda Cut-off parameter.
 *
 * Construct an exponentially cut off power law from
 * - a prefactor value (in units of ph/cm2/s/MeV)
 * - a spectral index,
 * - a pivot energy, and
 * - a cut-off parameter lambda.
 ***************************************************************************/
GModelSpectralExpPlaw2::GModelSpectralExpPlaw2(const double&  prefactor,
                                              const double&  index,
                                              const GEnergy& pivot,
                                              const double & lambda) :
                       GModelSpectral()
{
    // Initialise members
    init_members();

    // Set parameters
    m_norm.value(prefactor);
    m_index.value(index);
    m_pivot.value(pivot.MeV()); // Internally stored in MeV
    m_lambda.value(lambda);     // Internally stored in 1/MeV

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
GModelSpectralExpPlaw2::GModelSpectralExpPlaw2(const GXmlElement& xml) :
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
 * @param[in] model Exponentially cut off power law model.
 ***************************************************************************/
GModelSpectralExpPlaw2::GModelSpectralExpPlaw2(const GModelSpectralExpPlaw2& model) :
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
GModelSpectralExpPlaw2::~GModelSpectralExpPlaw2(void)
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
 * @param[in] model Exponentially cut off power law model.
 * @return Exponentially cut off power law model.
 ***************************************************************************/
GModelSpectralExpPlaw2& GModelSpectralExpPlaw2::operator=(const GModelSpectralExpPlaw2& model)
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
void GModelSpectralExpPlaw2::clear(void)
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
 * @brief Clone exponentially cut off power law model
 *
 * @return Pointer to deep copy of exponentially cut off power law model.
 ***************************************************************************/
GModelSpectralExpPlaw2* GModelSpectralExpPlaw2::clone(void) const
{
    // Clone exponentially cut off power law model
    return new GModelSpectralExpPlaw2(*this);
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
 *    \left( \frac{E}{\tt m\_pivot} \right)^{\tt m\_index}
 *    \exp \left( -{\tt m\_lamda}*E \right)
 * \f]
 *
 * where
 * - \f${\tt m\_norm}\f$ is the normalization or prefactor,
 * - \f${\tt m\_pivot}\f$ is the pivot energy,
 * - \f${\tt m\_index}\f$ is the spectral index, and
 * - \f${\tt m\_lambda}\f$ is the cut-off parameter.
 *
 * @todo The method expects that pivot!=0. Otherwise Inf or NaN
 *       may result. We should add a test that prevents using invalid
 *       values.
 ***************************************************************************/
double GModelSpectralExpPlaw2::eval(const GEnergy& srcEng,
                                   const GTime&   srcTime) const
{
    // Update the evaluation cache
    update_eval_cache(srcEng);

    // Compute function value
    double value = m_norm.value() * m_last_power;

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GModelSpectralExpPlaw::eval";
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
 *    \left( \frac{E}{\tt m\_pivot} \right)^{\tt m\_index}
 *    \exp \left( -{\tt m\_lamda}*E \right)
 * \f]
 *
 * where
 * - \f${\tt m\_norm}\f$ is the normalization or prefactor,
 * - \f${\tt m\_pivot}\f$ is the pivot energy,
 * - \f${\tt m\_index}\f$ is the spectral index, and
 * - \f${\tt m\_lambda}\f$ is the cut-off parameter.
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
 *      S_{\rm E}(E | t) \, \ln(E/{\tt m_pivot})
 * \f]
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta {\tt m\_lambda}} =
 *      S_{\rm E}(E | t) \, \left( -E \right)
 * \f]
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta {\tt m\_pivot}} =
 *      -S_{\rm E}(E | t) \,
 *      \left( \frac{{\tt m\_index}}{{\tt m\_pivot}} \right)
 * \f]
 *
 * @todo The method expects that pivot!=0. Otherwise Inf or NaN
 *       may result. We should add a test that prevents using invalid
 *       values.
 ***************************************************************************/
double GModelSpectralExpPlaw2::eval_gradients(const GEnergy& srcEng,
                                             const GTime&   srcTime)
{
    // Update the evaluation cache
    update_eval_cache(srcEng);

    // Compute function value
    double value = m_norm.value() * m_last_power;

    // Compute partial derivatives with respect to the parameter factor
    // values. The partial derivatives with respect to the parameter
    // values are obtained by division by the scale factor.
    double g_norm   = (m_norm.is_free())
                      ? m_norm.scale() * m_last_power : 0.0;
    double g_index  = (m_index.is_free())
                      ? value * m_index.scale() * std::log(m_last_e_norm) : 0.0;
    double g_lambda = (m_lambda.is_free())
                      ? -value * m_lambda.scale() * m_last_energy.MeV() : 0.0;
    double g_pivot  = (m_pivot.is_free())
                      ? -value * m_last_index / m_pivot.factor_value() : 0.0;

    // Set gradients
    m_norm.factor_gradient(g_norm);
    m_index.factor_gradient(g_index);
    m_lambda.factor_gradient(g_lambda);
    m_pivot.factor_gradient(g_pivot);

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GModelSpectralExpPlaw::eval_gradients";
        std::cout << "(srcEng=" << srcEng;
        std::cout << ", srcTime=" << srcTime << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", e_norm=" << m_last_e_norm;
        std::cout << ", e_lambda=" << m_last_e_lambda;
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
double GModelSpectralExpPlaw2::flux(const GEnergy& emin,
                                   const GEnergy& emax) const
{
    // Initialise flux
    double flux = 0.0;
    
    // Compute only if integration range is valid
    if (emin < emax) {

        // Setup integration kernel
        flux_kernel integrand(m_norm.value(),  m_index.value(),
                              m_pivot.value(), m_lambda.value());
        GIntegral integral(&integrand);

        // Get integration boundaries in MeV
        double e_min = emin.MeV();
        double e_max = emax.MeV();

        // Perform integration
        flux = integral.romberg(e_min, e_max);

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
double GModelSpectralExpPlaw2::eflux(const GEnergy& emin,
                                    const GEnergy& emax) const
{
    // Initialise flux
    double eflux = 0.0;
    
    // Compute only if integration range is valid
    if (emin < emax) {

        // Setup integration kernel
        eflux_kernel integrand(m_norm.value(),  m_index.value(),
                               m_pivot.value(), m_lambda.value());
        GIntegral integral(&integrand);

        // Get integration boundaries in MeV
        double e_min = emin.MeV();
        double e_max = emax.MeV();

        // Perform integration
        eflux = integral.romberg(e_min, e_max);

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
GEnergy GModelSpectralExpPlaw2::mc(const GEnergy& emin,
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
    double inv_ecut = m_lambda.value();

    // Use rejection method to draw a random energy. We first draw
    // analytically from a power law, and then compare the power law
    // at the drawn energy to the exponentially cut off function. This
    // gives an acceptance fraction, and we accept the energy only if
    // a uniform random number is <= the acceptance fraction.
    do {

        // Get uniform random number
        double u = ran.uniform();

        // Case A: Index is not -1
        if (index() != -1.0) {
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
        acceptance_fraction = std::exp(-eng * inv_ecut);

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
 *     <spectrum type="ExpCutoff2">
 *       <parameter name="Prefactor" scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Index"     scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="lambda"    scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Scale"     scale=".." value=".." min=".." max=".." free=".."/>
 *     </spectrum>
 ***************************************************************************/
void GModelSpectralExpPlaw2::read(const GXmlElement& xml)
{
    // Verify that XML element has exactly 4 parameters
    if (xml.elements() != 4 || xml.elements("parameter") != 4) {
        throw GException::model_invalid_parnum(G_READ, xml,
              "Power law model requires exactly 4 parameters.");
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

        // Handle index
        else if (par->attribute("name") == "Index") {
            m_index.read(*par);
            npar[1]++;
        }

        // Handle cutoff parameter
        else if (par->attribute("name") == "lambda") {
            m_lambda.read(*par);
            npar[2]++;
        }

        // Handle pivot energy
        else if (par->attribute("name") == "Scale") {
            m_pivot.read(*par);
            npar[3]++;
        }

    } // endfor: looped over all parameters

    // Verify that all parameters were found
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1 || npar[3] != 1) {
        throw GException::model_invalid_parnames(G_READ, xml,
              "Require \"Prefactor\", \"Index\", \"lambda\" and \"Scale\""
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
 *     <spectrum type="ExpCutoff2">
 *       <parameter name="Prefactor" scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Index"     scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="lambda"    scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Scale"     scale=".." value=".." min=".." max=".." free=".."/>
 *     </spectrum>
 ***************************************************************************/
void GModelSpectralExpPlaw2::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "") {
        xml.attribute("type", "ExpCutoff2");
    }

    // Verify model type
    if (xml.attribute("type") != "ExpCutoff2") {
        throw GException::model_invalid_spectral(G_WRITE, xml.attribute("type"),
              "Spectral model is not of type \"ExpCutoff2\".");
    }

    // If XML element has 0 nodes then append 4 parameter nodes
    if (xml.elements() == 0) {
        xml.append(GXmlElement("parameter name=\"Prefactor\""));
        xml.append(GXmlElement("parameter name=\"Index\""));
        xml.append(GXmlElement("parameter name=\"lambda\""));
        xml.append(GXmlElement("parameter name=\"Scale\""));
    }

    // Verify that XML element has exactly 4 parameters
    if (xml.elements() != 4 || xml.elements("parameter") != 4) {
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Power law model requires exactly 4 parameters.");
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

        // Handle index
        else if (par->attribute("name") == "Index") {
            npar[1]++;
            m_index.write(*par);
        }

        // Handle cutoff parameter
        else if (par->attribute("name") == "lambda") {
            npar[2]++;
            m_lambda.write(*par);
        }

        // Handle pivot energy
        else if (par->attribute("name") == "Scale") {
            m_pivot.write(*par);
            npar[3]++;
        }

    } // endfor: looped over all parameters

    // Check of all required parameters are present
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1 || npar[3] != 1) {
        throw GException::model_invalid_parnames(G_WRITE, xml,
              "Require \"Prefactor\", \"Index\", \"lambda\" and \"Scale\""
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
std::string GModelSpectralExpPlaw2::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpectralExpPlaw2 ===");

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
void GModelSpectralExpPlaw2::init_members(void)
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
    m_index.clear();
    m_index.name("Index");
    m_index.scale(1.0);
    m_index.value(-2.0);        // default: -2.0
    m_index.range(-10.0,+10.0); // range:   [-10,+10]
    m_index.free();
    m_index.gradient(0.0);
    m_index.has_grad(true);

    // Initialise cut off parameter
    m_lambda.clear();
    m_lambda.name("lambda");
    m_lambda.unit("1/MeV");
    m_lambda.scale(1.0);
    m_lambda.value(0.001);      // default: 0.001
    m_lambda.min(0.0);          // min:     0.0
    m_lambda.max(10.0);	    	// max:     10.0
    m_lambda.free();
    m_lambda.gradient(0.0);
    m_lambda.has_grad(true);

    // Initialise pivot energy
    m_pivot.clear();
    m_pivot.name("PivotEnergy");
    m_pivot.unit("MeV");
    m_pivot.scale(1.0);
    m_pivot.value(100.0);
    m_pivot.fix();
    m_pivot.gradient(0.0);
    m_pivot.has_grad(true);

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);
    m_pars.push_back(&m_index);
    m_pars.push_back(&m_lambda);
    m_pars.push_back(&m_pivot);

    // Initialise eval cache
    m_last_energy.clear();
    m_last_index    = 1.0e30;
    m_last_lambda   = 1.0e30;
    m_last_pivot    = 1.0e30;
    m_last_e_norm   = 0.0;
    m_last_e_lambda = 0.0;
    m_last_power    = 0.0;

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
void GModelSpectralExpPlaw2::copy_members(const GModelSpectralExpPlaw2& model)
{
    // Copy members
    m_norm   = model.m_norm;
    m_index  = model.m_index;
    m_lambda = model.m_lambda;
    m_pivot  = model.m_pivot;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);
    m_pars.push_back(&m_index);
    m_pars.push_back(&m_lambda);
    m_pars.push_back(&m_pivot);

    // Copy eval cache
    m_last_energy   = model.m_last_energy;
    m_last_index    = model.m_last_index;
    m_last_lambda   = model.m_last_lambda;
    m_last_pivot    = model.m_last_pivot;
    m_last_e_norm   = model.m_last_e_norm;
    m_last_e_lambda = model.m_last_e_lambda;
    m_last_power    = model.m_last_power;

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
void GModelSpectralExpPlaw2::free_members(void)
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
void GModelSpectralExpPlaw2::update_eval_cache(const GEnergy& energy) const
{
    // Get parameter values (takes 3 multiplications which are difficult
    // to avoid)
    double index  = m_index.value();
    double lambda = m_lambda.value();
    double pivot  = m_pivot.value();
    
    // If the energy or one of the parameters index, lambda or pivot
    // energy has changed then recompute the cache
    if ((m_last_energy != energy) ||
        (m_last_index  != index)  ||
        (m_last_lambda != lambda) ||
        (m_last_pivot  != pivot)) {

        // Store actual energy and parameter values
        m_last_energy = energy;
        m_last_index  = index;
        m_last_lambda = lambda;
        m_last_pivot  = pivot;

        // Compute and store value
        double eng      = energy.MeV();
        m_last_e_norm   = eng / m_last_pivot;
        m_last_e_lambda = eng * m_last_lambda;
        m_last_power    = std::pow(m_last_e_norm, m_last_index) *
                          std::exp(-m_last_e_lambda);
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
void GModelSpectralExpPlaw2::update_mc_cache(const GEnergy& emin,
                                             const GEnergy& emax) const

{
    // Case A: Index is not -1
    if (index() != -1.0) {

        // Change in energy boundaries?
        if (emin.MeV() != m_mc_emin || emax.MeV() != m_mc_emax) {
            m_mc_emin       = emin.MeV();
            m_mc_emax       = emax.MeV();
            m_mc_exponent   = index() + 1.0;
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
double GModelSpectralExpPlaw2::flux_kernel::eval(const double& energy)
{
    // Evaluate function value
    double e_norm = energy * m_inv_pivot;
    double e_cut  = energy * m_lambda;
    double power  = std::pow(e_norm, m_index) * std::exp(-e_cut);
    double value  = m_norm * power;

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Kernel for energy flux integration
 *
 * @param[in] energy Energy (MeV).
 ***************************************************************************/
double GModelSpectralExpPlaw2::eflux_kernel::eval(const double& energy)
{
    // Evaluate function value
    double e_norm = energy * m_inv_pivot;
    double e_cut  = energy * m_lambda;
    double power  = std::pow(e_norm, m_index) * std::exp(-e_cut);
    double value  = m_norm * power * energy;

    // Return value
    return value;
}
