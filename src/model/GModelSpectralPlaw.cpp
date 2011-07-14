/***************************************************************************
 *        GModelSpectralPlaw.cpp  -  Spectral power law model class        *
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
 * @file GModelSpectralPlaw.cpp
 * @brief Power law spectral model class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelSpectralPlaw.hpp"
#include "GModelSpectralRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpectralPlaw     g_spectral_plaw_seed;
const GModelSpectralRegistry g_spectral_plaw_registry(&g_spectral_plaw_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ                       "GModelSpectralPlaw::read(GXmlElement&)"
#define G_WRITE                     "GModelSpectralPlaw::write(GXmlElement&)"

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
GModelSpectralPlaw::GModelSpectralPlaw(void) : GModelSpectral()
{
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] norm Power law normalization.
 * @param[in] index Power law index.
 *
 * Construct a spectral power law from a normalization value and a spectral
 * index.
 ***************************************************************************/
GModelSpectralPlaw::GModelSpectralPlaw(const double& norm,
                                       const double& index) : GModelSpectral()
{
    // Initialise members
    init_members();

    // Set parameters
    m_norm.real_value(norm);
    m_index.real_value(index);

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element containing position information.
 *
 * Creates instance of a power law spectral model by extracting information
 * from an XML element. See GModelSpectralPlaw::read() for more information
 * about the expected structure of the XML element.
 ***************************************************************************/
GModelSpectralPlaw::GModelSpectralPlaw(const GXmlElement& xml) : GModelSpectral()
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
GModelSpectralPlaw::GModelSpectralPlaw(const GModelSpectralPlaw& model) :
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
GModelSpectralPlaw::~GModelSpectralPlaw(void)
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
 ***************************************************************************/
GModelSpectralPlaw& GModelSpectralPlaw::operator=(const GModelSpectralPlaw& model)
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
void GModelSpectralPlaw::clear(void)
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
GModelSpectralPlaw* GModelSpectralPlaw::clone(void) const
{
    return new GModelSpectralPlaw(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] srcEng True energy of photon.
 *
 * The power law function is defined as
 * \f[I(E)=norm (E/pivot)^{index}\f]
 * where
 * \f$norm\f$ is the normalization or prefactor,
 * \f$pivot\f$ is the pivot energy, and
 * \f$index\f$ is the spectral index.
 *
 * @todo For the moment the pivot energy is fixed to units of MeV. This may
 * not be ideal and should eventually be improved in the futur.
 * Furthermore, the method expects that energy!=0. Otherwise Inf or NaN
 * may result.
 ***************************************************************************/
double GModelSpectralPlaw::eval(const GEnergy& srcEng) const
{
    // Compute function value
    double energy = srcEng.MeV() / pivot();
    double power  = std::pow(energy, index());
    double value  = norm() * power;

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (isnan(value) || isinfinite(value)) {
        std::cout << "*** ERROR: GModelSpectralPlaw::eval";
        std::cout << "(srcEng=" << srcEng << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", energy=" << energy;
        std::cout << ", index=" << index();
        std::cout << ", pivot=" << pivot();
        std::cout << ", power=" << power;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] srcEng True energy of photon.
 *
 * The power law function is defined as
 * \f[I(E)=norm (E/pivot)^{index}\f]
 * where
 * \f$norm=n_s n_v\f$ is the normalization or prefactor,
 * \f$pivot=p_s p_v\f$ is the pivot energy, and
 * \f$index=i_s i_v\f$ is the spectral index.
 * Note that each parameter is factorised into a scaling factor and a value
 * and that the method is expected to return the gradient with respect to
 * the parameter value (i.e. n_v, p_v, and i_v in this case).
 *
 * The partial derivatives of the parameter values are given by
 * \f[dI/dn_v=n_s (E/pivot)^{index}\f]
 * \f[dI/dp_v=norm (E/pivot)^{index} (-index) / p_v\f]
 * \f[dI/di_v=norm (E/pivot)^{index} i_s \ln(E/pivot)\f]
 *
 * @todo For the moment the pivot energy is fixed to units of MeV. This may
 * not be ideal and should eventually be improved in the futur.
 * Furthermore, the method expects that energy!=0. Otherwise Inf or NaN
 * may result.
 ***************************************************************************/
double GModelSpectralPlaw::eval_gradients(const GEnergy& srcEng) const
{
    // Compute function value
    double energy = srcEng.MeV() / pivot();
    double power  = std::pow(energy, index());
    double value  = norm() * power;

    // Compute partial derivatives of the parameter values
    double g_norm  = (m_norm.isfree())  ? m_norm.scale() * power : 0.0;
    double g_index = (m_index.isfree()) ? value * m_index.scale() * std::log(energy) : 0.0;
    double g_pivot = (m_pivot.isfree()) ? -value * index() / m_pivot.value() : 0.0;

    // Set gradients (circumvent const correctness)
    ((GModelSpectralPlaw*)this)->m_norm.gradient(g_norm);
    ((GModelSpectralPlaw*)this)->m_index.gradient(g_index);
    ((GModelSpectralPlaw*)this)->m_pivot.gradient(g_pivot);

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (isnan(value) || isinfinite(value)) {
        std::cout << "*** ERROR: GModelSpectralPlaw::eval_gradients";
        std::cout << "(srcEng=" << srcEng << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", energy=" << energy;
        std::cout << ", index=" << index();
        std::cout << ", pivot=" << pivot();
        std::cout << ", power=" << power;
        std::cout << ", g_norm=" << g_norm;
        std::cout << ", g_index=" << g_index;
        std::cout << ", g_pivot=" << g_pivot;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Returns model flux between [emin, emax] (units: ph/cm2/s)
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
 ***************************************************************************/
double GModelSpectralPlaw::flux(const GEnergy& emin, const GEnergy& emax) const
{
    // Compute flux
    double flux = norm() * pow(pivot(), -index());
    if (index() != -1.0) {
        double exponent = index() + 1.0;
        flux *= (std::pow(emax.MeV(), exponent)-std::pow(emin.MeV(), exponent)) /
                exponent;
    }
    else
        flux *= (std::log(emax.MeV()) - std::log(emin.MeV()));

    // Return flux
    return flux;
}


/***********************************************************************//**
 * @brief Returns Monte Carlo energy between [emin, emax]
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 * @param[in] ran Random number generator.
 *
 * Returns Monte Carlo energy by randomly drawing from a power law.
 ***************************************************************************/
GEnergy GModelSpectralPlaw::mc(const GEnergy& emin, const GEnergy& emax,
                               GRan& ran) const
{
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
 * @brief Autoscale normalization
 *
 * Based on the actual value of the m_norm parameter, set the scale of m_norm
 * so that the value will be 1. If minimum and/or maximum value boundaries
 * exist, the boundaries are also modified accordingly.
 ***************************************************************************/
void GModelSpectralPlaw::autoscale(void)
{
    // Autoscale normalization to a value of 1.0
    if (m_norm.value() != 0.0) {

        // Get inverse scaling factor
        double invscale = 1.0 / m_norm.value();

        // Set values, error, min and max
        m_norm.value(m_norm.value() * invscale);
        m_norm.error(m_norm.error() * invscale);
        if (m_norm.hasmin())
            m_norm.min(m_norm.min() * invscale);
        if (m_norm.hasmax())
            m_norm.max(m_norm.max() * invscale);

        // Set scale
        m_norm.scale(1.0 / invscale);

    }

    // Return
    return;
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
 * Read the spectral power law information from an XML element. The XML
 * element is required to have 3 parameters with names "Prefactor", "Index",
 * and "Scale".
 *
 * @todo Add parameter validity check
 ***************************************************************************/
void GModelSpectralPlaw::read(const GXmlElement& xml)
{
    // Verify that XML element has exactly 3 parameters
    if (xml.elements() != 3 || xml.elements("parameter") != 3)
        throw GException::model_invalid_parnum(G_READ, xml,
              "Power law model requires exactly 3 parameters.");

    // Extract model parameters
    int npar[] = {0, 0, 0};
    for (int i = 0; i < 3; ++i) {

        // Get parameter element
        GXmlElement* par = (GXmlElement*)xml.element("parameter", i);

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

        // Handle pivot energy
        else if (par->attribute("name") == "Scale") {
            m_pivot.read(*par);
            npar[2]++;
        }

    } // endfor: looped over all parameters

    // Verify that all parameters were found
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1)
        throw GException::model_invalid_parnames(G_READ, xml,
              "Require \"Prefactor\", \"Index\" and \"Scale\" parameters.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * @exception GException::model_invalid_spectral
 *            Existing XML element is not of type "PowerLaw"
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters or nodes found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Write the spectral power law information into an XML element. The XML
 * element has to be of type "PowerLaw" and will have 3 parameter leafs
 * named "Prefactor", "Index", and "Scale".
 ***************************************************************************/
void GModelSpectralPlaw::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "")
        xml.attribute("type", "PowerLaw");

    // Verify model type
    if (xml.attribute("type") != "PowerLaw")
        throw GException::model_invalid_spectral(G_WRITE, xml.attribute("type"),
              "Spectral model is not of type \"PowerLaw\".");

    // If XML element has 0 nodes then append 3 parameter nodes
    if (xml.elements() == 0) {
        xml.append(new GXmlElement("parameter name=\"Prefactor\""));
        xml.append(new GXmlElement("parameter name=\"Index\""));
        xml.append(new GXmlElement("parameter name=\"Scale\""));
    }

    // Verify that XML element has exactly 3 parameters
    if (xml.elements() != 3 || xml.elements("parameter") != 3)
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Power law model requires exactly 3 parameters.");

    // Set or update model parameter attributes
    int npar[] = {0, 0, 0};
    for (int i = 0; i < 3; ++i) {

        // Get parameter element
        GXmlElement* par = (GXmlElement*)xml.element("parameter", i);

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

        // Handle pivot energy
        else if (par->attribute("name") == "Scale") {
            m_pivot.write(*par);
            npar[2]++;
        }

    } // endfor: looped over all parameters

    // Check of all required parameters are present
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1)
        throw GException::model_invalid_parnames(G_WRITE, xml,
              "Require \"Prefactor\", \"Index\" and \"Scale\" parameters.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print powerlaw information
 ***************************************************************************/
std::string GModelSpectralPlaw::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GModelSpectralPlaw ===\n");
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
void GModelSpectralPlaw::init_members(void)
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

    // Initialise powerlaw index
    m_index.clear();
    m_index.name("Index");
    m_index.scale(1.0);
    m_index.value(-2.0);        // default: -2.0
    m_index.range(-10.0,+10.0); // range:   [-10,+10]
    m_index.free();
    m_index.gradient(0.0);
    m_index.hasgrad(true);

    // Initialise pivot energy
    m_pivot.clear();
    m_pivot.name("PivotEnergy");
    m_pivot.unit("MeV");
    m_pivot.scale(1.0);
    m_pivot.value(100.0);       // default: 100
    m_pivot.fix();
    m_pivot.gradient(0.0);
    m_pivot.hasgrad(true);

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);
    m_pars.push_back(&m_index);
    m_pars.push_back(&m_pivot);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model GModelSpectralPlaw members which should be copied.
 ***************************************************************************/
void GModelSpectralPlaw::copy_members(const GModelSpectralPlaw& model)
{
    // Copy members
    m_norm  = model.m_norm;
    m_index = model.m_index;
    m_pivot = model.m_pivot;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);
    m_pars.push_back(&m_index);
    m_pars.push_back(&m_pivot);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralPlaw::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
