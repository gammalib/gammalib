/***************************************************************************
 *         GModelSpectralPlaw2.cpp - Spectral power law model class        *
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
 * @file GModelSpectralPlaw2.cpp
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
#include "GModelSpectralPlaw2.hpp"
#include "GModelSpectralRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpectralPlaw2    g_spectral_plaw2_seed;
const GModelSpectralRegistry g_spectral_plaw2_registry(&g_spectral_plaw2_seed);

/* __ Method name definitions ____________________________________________ */
#define G_FLUX                "GModelSpectralPlaw2::flux(GEnergy&, GEnergy&)"
#define G_EFLUX              "GModelSpectralPlaw2::eflux(GEnergy&, GEnergy&)"
#define G_MC     "GModelSpectralPlaw2::mc(GEnergy&, GEnergy&, GTime&, GRan&)"
#define G_READ                      "GModelSpectralPlaw2::read(GXmlElement&)"
#define G_WRITE                    "GModelSpectralPlaw2::write(GXmlElement&)"

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
GModelSpectralPlaw2::GModelSpectralPlaw2(void) : GModelSpectral()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] integral Integral flux (ph/cm2/s).
 * @param[in] index Power law index.
 * @param[in] emin Minimum energy.
 * @param[in] emax Maximum energy.
 *
 * Construct a spectral power law from the
 * - integral flux (in ph/cm2/s),
 * - spectral index,
 * - minimum energy and
 * - maximum energy.
 ***************************************************************************/
GModelSpectralPlaw2::GModelSpectralPlaw2(const double&  integral,
                                         const double&  index,
                                         const GEnergy& emin,
                                         const GEnergy& emax) :
                     GModelSpectral()
{
    // Initialise members
    init_members();

    // Set parameters
    m_integral.value(integral);
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
GModelSpectralPlaw2::GModelSpectralPlaw2(const GXmlElement& xml) :
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
GModelSpectralPlaw2::GModelSpectralPlaw2(const GModelSpectralPlaw2& model) :
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
GModelSpectralPlaw2::~GModelSpectralPlaw2(void)
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
GModelSpectralPlaw2& GModelSpectralPlaw2::operator=(const GModelSpectralPlaw2& model)
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
void GModelSpectralPlaw2::clear(void)
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
GModelSpectralPlaw2* GModelSpectralPlaw2::clone(void) const
{
    // Clone power law model
    return new GModelSpectralPlaw2(*this);
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
 *    S_{\rm E}(E | t) = {\tt m\_integral}
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
 *    \frac{{\tt m\_integral}}
 *         {\log {\tt e\_max} - \log {\tt e\_min}}
 *    E^{\tt m\_index}
 * \f]
 *
 * for \f${\tt m\_index} = -1\f$, where
 * - \f${\tt e\_min}\f$ is the minimum energy of an interval,
 * - \f${\tt e\_max}\f$ is the maximum energy of an interval,
 * - \f${\tt m\_integral}\f$ is the integral flux between 
 *   \f${\tt e\_min}\f$ and \f${\tt e\_max}\f$, and
 * - \f${\tt m\_index}\f$ is the spectral index.
 ***************************************************************************/
double GModelSpectralPlaw2::eval(const GEnergy& srcEng,
                                 const GTime&   srcTime) const
{
    // Update precomputed values
    update(srcEng);

    // Compute function value
    double value = m_integral.value() * m_norm * m_power;

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (isnotanumber(value) || isinfinite(value)) {
        std::cout << "*** ERROR: GModelSpectralPlaw2::eval";
        std::cout << "(srcEng=" << srcEng;
        std::cout << ", srcTime=" << srcTime << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", integral=" << integral();
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
 *    S_{\rm E}(E | t) = {\tt m\_integral}
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
 *    \frac{{\tt m\_integral}}
 *         {\log {\tt e\_max} - \log {\tt e\_min}}
 *    E^{\tt m\_index}
 * \f]
 *
 * for \f${\tt m\_index} = -1\f$, where
 * - \f${\tt e\_min}\f$ is the minimum energy of an interval,
 * - \f${\tt e\_max}\f$ is the maximum energy of an interval,
 * - \f${\tt m\_integral}\f$ is the integral flux between 
 *   \f${\tt e\_min}\f$ and \f${\tt e\_max}\f$, and
 * - \f${\tt m\_index}\f$ is the spectral index.
 *
 * The method also evaluates the partial derivatives of the model with
 * respect to the parameters using
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta {\tt m\_integral}} =
 *      \frac{S_{\rm E}(E | t)}{{\tt m\_integral}}
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
double GModelSpectralPlaw2::eval_gradients(const GEnergy& srcEng,
                                           const GTime&   srcTime)
{
    // Initialise gradients
    double g_integral = 0.0;
    double g_index    = 0.0;
    
    // Update precomputed values
    update(srcEng);

    // Compute function value
    double value = m_integral.value() * m_norm * m_power;

    // Integral flux gradient
    if (m_integral.isfree()) {
         g_integral = value / m_integral.factor_value();
    }

    // Index gradient
    if (m_index.isfree()) {
        g_index = value * (m_g_norm + ln10*srcEng.log10MeV()) * m_index.scale();
    }

    // Set gradients
    m_integral.factor_gradient(g_integral);
    m_index.factor_gradient(g_index);

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (isnotanumber(value) || isinfinite(value)) {
        std::cout << "*** ERROR: GModelSpectralPlaw2::eval_gradients";
        std::cout << "(srcEng=" << srcEng;
        std::cout << ", srcTime=" << srcTime << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", integral=" << integral();
        std::cout << ", m_norm=" << m_norm;
        std::cout << ", m_power=" << m_power;
        std::cout << ", g_integral=" << g_integral;
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
double GModelSpectralPlaw2::flux(const GEnergy& emin,
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
            flux                = m_integral.value() * factor;
        }

        // Case B: Index is -1
        else {
            double log_emin     = std::log(emin.MeV());
            double log_emax     = std::log(emax.MeV());
            double log_ref_emin = std::log(this->emin().MeV());
            double log_ref_emax = std::log(this->emax().MeV());
            double factor       = (log_emax - log_emin) /
                                  (log_ref_emax - log_ref_emin);
            flux                = m_integral.value() * factor;
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
double GModelSpectralPlaw2::eflux(const GEnergy& emin,
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
            norm = m_integral.value() * gamma / (pow_ref_emax - pow_ref_emin);
        }
        else {
            double log_ref_emin = std::log(this->emin().MeV());
            double log_ref_emax = std::log(this->emax().MeV());
            norm = m_integral.value() / (log_ref_emax - log_ref_emin);
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
        eflux *= MeV2erg;

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
GEnergy GModelSpectralPlaw2::mc(const GEnergy& emin,
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
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Reads the spectral information from an XML element. The format of the XML
 * elements is
 *
 *     <spectrum type="PowerLaw2">
 *       <parameter name="Integral"   scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Index"      scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="LowerLimit" scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="UpperLimit" scale=".." value=".." min=".." max=".." free=".."/>
 *     </spectrum>
 *
 * @todo Add parameter validity check
 ***************************************************************************/
void GModelSpectralPlaw2::read(const GXmlElement& xml)
{
    // Verify that XML element has exactly 4 parameters
    if (xml.elements() != 4 || xml.elements("parameter") != 4) {
        throw GException::model_invalid_parnum(G_READ, xml,
              "Power law 2 spectral model requires exactly 4 parameters.");
    }

    // Extract model parameters
    int npar[] = {0, 0, 0, 0};
    for (int i = 0; i < 4; ++i) {

        // Get parameter element
        const GXmlElement* par = xml.element("parameter", i);

        // Handle prefactor
        if (par->attribute("name") == "Integral") {
            m_integral.read(*par);
            npar[0]++;
        }

        // Handle index
        else if (par->attribute("name") == "Index") {
            m_index.read(*par);
            npar[1]++;
        }

        // Handle lower limit
        else if (par->attribute("name") == "LowerLimit") {
            m_emin.read(*par);
            npar[2]++;
        }

        // Handle upper limit
        else if (par->attribute("name") == "UpperLimit") {
            m_emax.read(*par);
            npar[3]++;
        }

    } // endfor: looped over all parameters

    // Verify that all parameters were found
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1 || npar[3] != 1) {
        throw GException::model_invalid_parnames(G_READ, xml,
              "Power law 2 spectral model requires \"Integral\", \"Index\","
              " \"LowerLimit\" and \"UpperLimit\" parameters.");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::model_invalid_spectral
 *            Existing XML element is not of type "PowerLaw2"
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters or nodes found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Writes the spectral information into an XML element. The format of the XML
 * element is
 *
 *     <spectrum type="PowerLaw2">
 *       <parameter name="Integral"   scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Index"      scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="LowerLimit" scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="UpperLimit" scale=".." value=".." min=".." max=".." free=".."/>
 *     </spectrum>
 ***************************************************************************/
void GModelSpectralPlaw2::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "") {
        xml.attribute("type", "PowerLaw2");
    }

    // Verify model type
    if (xml.attribute("type") != "PowerLaw2") {
        throw GException::model_invalid_spectral(G_WRITE, xml.attribute("type"),
              "Spectral model is not of type \"PowerLaw2\".");
    }

    // If XML element has 0 nodes then append 4 parameter nodes
    if (xml.elements() == 0) {
        xml.append(GXmlElement("parameter name=\"Integral\""));
        xml.append(GXmlElement("parameter name=\"Index\""));
        xml.append(GXmlElement("parameter name=\"LowerLimit\""));
        xml.append(GXmlElement("parameter name=\"UpperLimit\""));
    }

    // Verify that XML element has exactly 4 parameters
    if (xml.elements() != 4 || xml.elements("parameter") != 4) {
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Power law 2 spectral model requires exactly 4 parameters.");
    }

    // Set or update model parameter attributes
    int npar[] = {0, 0, 0, 0};
    for (int i = 0; i < 4; ++i) {

        // Get parameter element
        GXmlElement* par = xml.element("parameter", i);

        // Handle prefactor
        if (par->attribute("name") == "Integral") {
            npar[0]++;
            m_integral.write(*par);
        }

        // Handle index
        else if (par->attribute("name") == "Index") {
            npar[1]++;
            m_index.write(*par);
        }

        // Handle lower limit
        else if (par->attribute("name") == "LowerLimit") {
            m_emin.write(*par);
            npar[2]++;
        }

        // Handle lower limit
        else if (par->attribute("name") == "UpperLimit") {
            m_emax.write(*par);
            npar[3]++;
        }

    } // endfor: looped over all parameters

    // Check of all required parameters are present
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1 || npar[3] != 1) {
        throw GException::model_invalid_parnames(G_WRITE, xml,
              "Power law 2 spectral model requires \"Integral\", \"Index\","
              " \"LowerLimit\" and \"UpperLimit\" parameters.");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print power law information
 *
 * @return String containing power law information.
 ***************************************************************************/
std::string GModelSpectralPlaw2::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GModelSpectralPlaw2 ===\n");

    // Append information
    result.append(parformat("Number of parameters")+str(size()));
    for (int i = 0; i < size(); ++i) {
        result.append("\n"+m_pars[i]->print());
    }

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
void GModelSpectralPlaw2::init_members(void)
{
    // Initialise integral flux
    m_integral.clear();
    m_integral.name("Integral");
    m_integral.unit("ph/cm2/s");
    m_integral.scale(1.0);
    m_integral.value(1.0);       // default: 1.0
    m_integral.range(0.0, 10.0); // range:   [0,10]
    m_integral.free();
    m_integral.gradient(0.0);
    m_integral.hasgrad(true);

    // Initialise powerlaw index
    m_index.clear();
    m_index.name("Index");
    m_index.scale(1.0);
    m_index.value(-2.0);        // default: -2.0
    m_index.range(-10.0,+10.0); // range:   [-10,+10]
    m_index.free();
    m_index.gradient(0.0);
    m_index.hasgrad(true);

    // Initialise lower limit
    m_emin.clear();
    m_emin.name("LowerLimit");
    m_emin.unit("MeV");
    m_emin.scale(1.0);
    m_emin.value(100.0);         // default: 100
    m_emin.range(0.001, 1.0e15); // range:   [0.001, 1e15]
    m_emin.fix();
    m_emin.gradient(0.0);
    m_emin.hasgrad(false);

    // Initialise upper limit
    m_emax.clear();
    m_emax.name("UpperLimit");
    m_emax.unit("MeV");
    m_emax.scale(1.0);
    m_emax.value(500000.0);      // default: 500000
    m_emin.range(0.001, 1.0e15); // range:   [0.001, 1e15]
    m_emax.fix();
    m_emax.gradient(0.0);
    m_emax.hasgrad(false);

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_integral);
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
    m_last_integral   = 0.0;
    m_last_index      = 1000.0;
    m_last_emin.MeV(0.0);
    m_last_emax.MeV(0.0);
    m_last_energy.MeV(0.0);
    m_last_value      = 0.0;
    m_last_g_integral = 0.0;
    m_last_g_index    = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Spectral power law model.
 ***************************************************************************/
void GModelSpectralPlaw2::copy_members(const GModelSpectralPlaw2& model)
{
    // Copy members
    m_integral = model.m_integral;
    m_index    = model.m_index;
    m_emin     = model.m_emin;
    m_emax     = model.m_emax;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_integral);
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
    m_last_integral   = model.m_last_integral;
    m_last_index      = model.m_last_index;
    m_last_emin       = model.m_last_emin;
    m_last_emax       = model.m_last_emax;
    m_last_energy     = model.m_last_energy;
    m_last_value      = model.m_last_value;
    m_last_g_integral = model.m_last_g_integral;
    m_last_g_index    = model.m_last_g_index;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralPlaw2::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update precomputed values
 *
 * @param[in] srcEng Source energy
 ***************************************************************************/
void GModelSpectralPlaw2::update(const GEnergy& srcEng) const
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
