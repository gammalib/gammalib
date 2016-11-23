/***************************************************************************
 *    GModelSpectralLogParabola.cpp - Log parabola spectral model class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2016 by Michael Mayer                               *
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
 * @file GModelSpectralLogParabola.cpp
 * @brief Log parabola spectral model class definition
 * @author Michael Mayer 
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
#include "GModelSpectralLogParabola.hpp"
#include "GModelSpectralRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpectralLogParabola g_spectral_logp_seed1("LogParabola",
                                                      "Prefactor",
                                                      "Index",
                                                      "PivotEnergy",
                                                      "Curvature");
const GModelSpectralRegistry    g_spectral_logp_registry1(&g_spectral_logp_seed1);
#if defined(G_LEGACY_XML_FORMAT)
const GModelSpectralLogParabola g_spectral_logp_seed2("LogParabola",
                                                      "norm",
                                                      "alpha",
                                                      "Eb",
                                                      "beta");
const GModelSpectralLogParabola g_spectral_logp_seed3("LogParabola",
                                                      "Prefactor",
                                                      "Index",
                                                      "Scale",
                                                      "Curvature");
const GModelSpectralRegistry    g_spectral_logp_registry2(&g_spectral_logp_seed2);
const GModelSpectralRegistry    g_spectral_logp_registry3(&g_spectral_logp_seed3);
#endif

/* __ Method name definitions ____________________________________________ */
#define G_FLUX          "GModelSpectralLogParabola::flux(GEnergy&, GEnergy&)"
#define G_EFLUX        "GModelSpectralLogParabola::eflux(GEnergy&, GEnergy&)"
#define G_MC      "GModelSpectralLogParabola::mc(GEnergy&, GEnergy&, GTime&,"\
                                                                    " GRan&)"
#define G_READ                "GModelSpectralLogParabola::read(GXmlElement&)"
#define G_WRITE              "GModelSpectralLogParabola::write(GXmlElement&)"

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
GModelSpectralLogParabola::GModelSpectralLogParabola(void) : GModelSpectral()
{
    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Model type and parameter name constructor
 *
 * @param[in] type Model type.
 * @param[in] prefactor Name of prefactor parameter.
 * @param[in] index Name of index parameter.
 * @param[in] pivot Name of pivot parameter.
 * @param[in] curvature Name of curvature parameter.
 ***************************************************************************/
GModelSpectralLogParabola::GModelSpectralLogParabola(const std::string& type,
                                                     const std::string& prefactor,
                                                     const std::string& index,
                                                     const std::string& pivot,
                                                     const std::string& curvature) :
                           GModelSpectral()
{
    // Initialise members
    init_members();

    // Set model type
    m_type = type;

    // Set parameter names
    m_norm.name(prefactor);
    m_index.name(index);
    m_pivot.name(pivot);
    m_curvature.name(curvature);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] prefactor Power law pre factor.
 * @param[in] index Power law index.
 * @param[in] pivot Pivot energy.
 * @param[in] curvature Curvature.
 *
 * Construct a LogParabola model from
 * - a pre factor (in ph/cm2/s/MeV),
 * - a spectral index,
 * - a pivot energy, and
 * - a curvature.
 ***************************************************************************/
GModelSpectralLogParabola::GModelSpectralLogParabola(const double&  prefactor,
                                                     const double&  index,
                                                     const GEnergy& pivot,
                                                     const double&  curvature) :
                           GModelSpectral()
{
    // Initialise members
    init_members();

    // Set parameters
    m_norm.value(prefactor);
    m_index.value(index);
    m_pivot.value(pivot.MeV());
    m_curvature.value(curvature);

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
 * Constructs log parabola spectral model by extracting information from an
 * XML element. See the read() method for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelSpectralLogParabola::GModelSpectralLogParabola(const GXmlElement& xml) : 
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
 * @param[in] model LogParabola model.
 ***************************************************************************/
GModelSpectralLogParabola::GModelSpectralLogParabola(const GModelSpectralLogParabola& model) :
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
GModelSpectralLogParabola::~GModelSpectralLogParabola(void)
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
 * @param[in] model LogParabola model.
 * @return LogParabola model.
 ***************************************************************************/
GModelSpectralLogParabola& GModelSpectralLogParabola::operator=(const GModelSpectralLogParabola& model)
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
 * @brief Clear log parabola model
 ***************************************************************************/
void GModelSpectralLogParabola::clear(void)
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
 * @brief Clone log parabola model
 *
 * @return Pointer to deep copy of log parabola spectrum.
 ***************************************************************************/
GModelSpectralLogParabola* GModelSpectralLogParabola::clone(void) const
{
    // Clone log parabola model
    return new GModelSpectralLogParabola(*this);
}


/***********************************************************************//**
 * @brief Evaluate model
 *
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @param[in] gradients Compute gradients?
 * @return Model value (ph/cm2/s/MeV).
 *
 * Computes
 *
 * \f[
 *    S_{\rm E}(E | t) = {\tt m\_norm}
 *    \left( \frac{E}{\tt m\_pivot} \right)^{{\tt m\_index} +
 *    {\tt m\_curvature} \, \ln \frac{E}{\tt m\_pivot}}
 * \f]
 *
 * where
 * - \f${\tt m\_norm}\f$ is the normalization or prefactor,
 * - \f${\tt m\_index}\f$ is the spectral index,
 * - \f${\tt m\_curvature}\f$ is the spectral curvature, and
 * - \f${\tt m\_pivot}\f$ is the pivot energy.
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
 *      S_{\rm E}(E | t) \, \ln(E/{\tt m_pivot})
 * \f]
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta {\tt m\_curvature}} =
 *      S_{\rm E}(E | t) \, (\ln(E/{\tt m_pivot})^2)
 * \f]
 *
 * \f[
 *    \frac{\delta S_{\rm E}(E | t)}{\delta {\tt m\_pivot}} =
 *      -S_{\rm E}(E | t) \,
 *      \left( \frac{2 {\tt m\_curvature} \ln(E/{\tt m_pivot}) + {\tt m\_index}}
 *                  {{\tt m\_pivot}} \right)
 * \f]
 *
 * @todo The method expects that energy!=0. Otherwise Inf or NaN may result.
 ***************************************************************************/
double GModelSpectralLogParabola::eval(const GEnergy& srcEng,
                                       const GTime&   srcTime,
                                       const bool&    gradients) const
{
    // Update the evaluation cache
    update_eval_cache(srcEng);

    // Compute function value
    double value = m_norm.value() * m_last_power;

    // Optionally compute gradients
    if (gradients) {

        // Compute partial derivatives of the parameter values
        double g_norm  = (m_norm.is_free())
                         ? m_norm.scale() * m_last_power : 0.0;
        double g_index = (m_index.is_free())
                         ? value * m_index.scale() * m_last_log_e_norm : 0.0;
        double g_curvature = (m_curvature.is_free())
                             ? value * m_curvature.scale() * m_last_log_e_norm *
                             m_last_log_e_norm : 0.0;
        double g_pivot = (m_pivot.is_free())
                         ? -value * (m_last_exponent + m_curvature.value() *
                         m_last_log_e_norm) / m_pivot.factor_value() : 0.0;

        // Set gradients
        m_norm.factor_gradient(g_norm);
        m_index.factor_gradient(g_index);
        m_curvature.factor_gradient(g_curvature);
        m_pivot.factor_gradient(g_pivot);

    } // endif: gradient computation was requested

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GModelSpectralLogParabola::eval";
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
double GModelSpectralLogParabola::flux(const GEnergy& emin,
                                       const GEnergy& emax) const
{
    // Initialise flux
    double flux = 0.0;
    
    // Compute only if integration range is valid
    if (emin < emax) {

        // Initialise function to integrate
        flux_kern kernel(prefactor(), index(), curvature(), pivot());

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
double GModelSpectralLogParabola::eflux(const GEnergy& emin,
                                        const GEnergy& emax) const
{
    // Initialise flux
    double eflux = 0.0;
    
    // Compute only if integration range is valid
    if (emin < emax) {
	
        // Initialise function to integrate
        eflux_kern kernel(prefactor(), index(), curvature(), pivot());

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
 * Returns Monte Carlo energy by randomly drawing from the spectral model.
 ***************************************************************************/
GEnergy GModelSpectralLogParabola::mc(const GEnergy& emin,
                                      const GEnergy& emax,
                                      const GTime&   time,
                                      GRan&          ran) const
{
    if (emin >= emax) {
	    throw GException::erange_invalid(G_MC, emin.MeV(), emax.MeV(),
	          "Minimum energy < maximum energy required.");
	}

	// Allocate energy
	GEnergy energy;

	// Update cache
	update_mc_cache(emin, emax, time);

	// Initialise energy
	double eng;

	// Initialse acceptance fraction
	double acceptance_fraction;

	// Use rejection method to draw a random energy. We first draw
	// analytically from a power law, and then compare the power law
	// at the drawn energy to the curved function. This
	// gives an acceptance fraction, and we accept the energy only if
	// a uniform random number is <= the acceptance fraction.
	do {

	    // Get uniform random number
	    double u = ran.uniform();

	    // Case A: Corresponding mc Plaw-Index is not -1
	    if (m_mc_exponent != 0.0) {
	        if (u > 0.0) {
	            eng = std::exp(std::log(u * m_mc_pow_ewidth + m_mc_pow_emin) /
	                           m_mc_exponent);
	        }
	        else {
	            eng = 0.0;
	        }
	    }

	    // Case B: Corresponding mc Plaw-Index is  -1
	    else {
	        eng = std::exp(u * m_mc_pow_ewidth + m_mc_pow_emin);
	    }

	    // Compute powerlaw at given energy
	    double e_norm = eng / m_pivot.value();
	    double plaw   = m_mc_norm * std::pow(e_norm, m_mc_exponent-1.0);

	    // Compute logparabola at given energy
	    double logparabola = prefactor() * 
                             std::pow(e_norm,index()+curvature()*std::log(e_norm));

	    // Compute acceptance fraction
	    acceptance_fraction = logparabola / plaw;

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
 * Read the spectral log parabola information from an XML element.
 ***************************************************************************/
void GModelSpectralLogParabola::read(const GXmlElement& xml)
{
    // Get parameter pointers
    const GXmlElement* norm      = gammalib::xml_get_par(G_READ, xml, m_norm.name());
    const GXmlElement* index     = gammalib::xml_get_par(G_READ, xml, m_index.name());
    const GXmlElement* curvature = gammalib::xml_get_par(G_READ, xml, m_curvature.name());
    const GXmlElement* pivot     = gammalib::xml_get_par(G_READ, xml, m_pivot.name());

    // Read parameters
    m_norm.read(*norm);
    m_index.read(*index);
    m_curvature.read(*curvature);
    m_pivot.read(*pivot);

    // In case that legacy parameters were used we need to negate the index
    // and curvature because they are differently defined
    if (gammalib::xml_has_par(xml, "alpha") &&
        gammalib::xml_has_par(xml, "beta")) {
        m_index.min(-m_index.max());
        m_index.max(-m_index.min());
        m_index.value(-m_index.value());
        m_curvature.min(-m_curvature.max());
        m_curvature.max(-m_curvature.min());
        m_curvature.value(-m_curvature.value());
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
 *            Existing XML element is not of type "LogParabola"
 *
 * Write the LogParabola model information into an XML element.
 ***************************************************************************/
void GModelSpectralLogParabola::write(GXmlElement& xml) const
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

    // Check if parameters are Fermi/LAT style and need to be negated
    GModelPar inx = m_index;
    GModelPar crv = m_curvature;
    if (inx.name() == "alpha") {
        inx.min(-inx.max());
        inx.max(-inx.min());
        inx.value(-inx.value());
    }
    if (crv.name() == "beta") {
        crv.min(-crv.max());
        crv.max(-crv.min());
        crv.value(-crv.value());
    }

    // Get XML parameters
    GXmlElement* norm      = gammalib::xml_need_par(G_WRITE, xml, m_norm.name());
    GXmlElement* index     = gammalib::xml_need_par(G_WRITE, xml, inx.name());
    GXmlElement* curvature = gammalib::xml_need_par(G_WRITE, xml, crv.name());
    GXmlElement* pivot     = gammalib::xml_need_par(G_WRITE, xml, m_pivot.name());

    // Write parameters
    m_norm.write(*norm);
    inx.write(*index);
    crv.write(*curvature);
    m_pivot.write(*pivot);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print LogParabola information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing LogParabola information.
 ***************************************************************************/
std::string GModelSpectralLogParabola::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpectralLogParabola ===");

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
void GModelSpectralLogParabola::init_members(void)
{
    // Initialise model type
    m_type = "LogParabola";

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

    // Initialise curvature
    m_curvature.clear();
    m_curvature.name("Curvature");
    m_curvature.scale(1.0);
    m_curvature.value(-0.1);        // default: -2.0
    m_curvature.range(-10.0,+10.0); // range:   [-10,+10]
    m_curvature.free();
    m_curvature.gradient(0.0);
    m_curvature.has_grad(true);

    // Initialise pivot energy
    m_pivot.clear();
    m_pivot.name("PivotEnergy");
    m_pivot.unit("MeV");
    m_pivot.scale(1.0);
    m_pivot.value(100.0);       // default: 100
    m_pivot.fix();
    m_pivot.gradient(0.0);
    m_pivot.has_grad(true);

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);
    m_pars.push_back(&m_index);
    m_pars.push_back(&m_curvature);
    m_pars.push_back(&m_pivot);

    // Initialise eval cache
    m_last_energy.clear();
    m_last_index      = 1.0e30;
    m_last_curvature  = 1.0e30;
    m_last_pivot      = 1.0e30;
    m_last_e_norm     = 0.0;
    m_last_log_e_norm = 0.0;
    m_last_exponent   = 0.0;
    m_last_power      = 0.0;

    // Initialise MC cache
    m_mc_emin       = 0.0;
    m_mc_emax       = 0.0;
    m_mc_exponent   = 0.0;
    m_mc_pow_emin   = 0.0;
    m_mc_pow_ewidth = 0.0;
    m_mc_norm       = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model GModelSpectralLogParabola members which should be copied.
 ***************************************************************************/
void GModelSpectralLogParabola::copy_members(const GModelSpectralLogParabola& model)
{
    // Copy members
    m_type      = model.m_type;
    m_norm      = model.m_norm;
    m_index     = model.m_index;
    m_curvature = model.m_curvature;
    m_pivot     = model.m_pivot;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);
    m_pars.push_back(&m_index);
    m_pars.push_back(&m_curvature);
    m_pars.push_back(&m_pivot);

    // Copy eval cache
    m_last_energy     = model.m_last_energy;
    m_last_index      = model.m_last_index;
    m_last_curvature  = model.m_last_curvature;
    m_last_pivot      = model.m_last_pivot;
    m_last_e_norm     = model.m_last_e_norm;
    m_last_log_e_norm = model.m_last_log_e_norm;
    m_last_exponent   = model.m_last_exponent;
    m_last_power      = model.m_last_power;

    // Copy MC cache
    m_mc_emin       = model.m_mc_emin;
    m_mc_emax       = model.m_mc_emax;
    m_mc_exponent   = model.m_mc_exponent;
    m_mc_pow_emin   = model.m_mc_pow_emin;
    m_mc_pow_ewidth = model.m_mc_pow_ewidth;
    m_mc_norm       = model.m_mc_norm;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralLogParabola::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update eval precomputation cache
 *
 * @param[in] energy Energy.
 *
 * Updates the precomputation cache for eval() method.
 ***************************************************************************/
void GModelSpectralLogParabola::update_eval_cache(const GEnergy& energy) const
{
    // Get parameter values (takes 3 multiplications which are difficult
    // to avoid)
    double index     = m_index.value();
    double curvature = m_curvature.value();
    double pivot     = m_pivot.value();
    
    // If the energy or one of the parameters index, curvature or pivot energy
    // has changed then recompute the cache
    if ((m_last_energy    != energy)    ||
        (m_last_index     != index)     ||
        (m_last_curvature != curvature) ||
        (m_last_pivot     != pivot)) {

        // Store actual energy and parameter values
        m_last_energy    = energy;
        m_last_index     = index;
        m_last_curvature = curvature;
        m_last_pivot     = pivot;

        // Compute and store value
        double eng        = energy.MeV();
        m_last_e_norm     = eng / m_last_pivot;
        m_last_log_e_norm = std::log(m_last_e_norm);
        m_last_exponent   = m_last_index + m_last_curvature * m_last_log_e_norm;
        m_last_power      = std::pow(m_last_e_norm, m_last_exponent);

    } // endif: recomputation was required

    // Return
    return;
}


/***********************************************************************//**
 * @brief Update Monte Carlo pre computation cache
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 * @param[in] time True photon arrival time.
 *
 * Updates the precomputation cache for Monte Carlo simulations.
 ***************************************************************************/
void GModelSpectralLogParabola::update_mc_cache(const GEnergy& emin,
                                                const GEnergy& emax,
                                                const GTime&   time) const

{
	// Only update if boundaries have changed
	if(emin.MeV() != m_mc_emin || emax.MeV() != m_mc_emax){

        // Store energy boundaries
		m_mc_emin = emin.MeV();
		m_mc_emax = emax.MeV();

		// Find a corresponding power law with the criterion
        // Plaw > LogParabola in the given interval
	    double index_pl;

		// Checking the sign of spectrum curvature
		if (curvature() < 0) {

			// Use the spectral index at the pivot energy of the LogParabola
			index_pl  = index();
			m_mc_norm = prefactor();
		}
		else {
			// Use a power law which connects the ends of the convex,
            // curved model

			// Plaw index defined by the slope of a straight line in the
            // log-log-plane
			index_pl = std::log(eval(emin,time)/eval(emax,time)) / 
                       std::log(emin.MeV()/emax.MeV());

			// Plaw norm defined such that Plaw = LogParabola at emin
			m_mc_norm = eval(emin,time) / std::pow(emin.MeV()/m_pivot.value(),index_pl);

		}

        // Set precomputation cache
		if (index_pl != -1.0) {
			m_mc_exponent   = index_pl + 1.0;
			m_mc_pow_emin   = std::pow(emin.MeV(), m_mc_exponent);
			m_mc_pow_ewidth = std::pow(emax.MeV(), m_mc_exponent) - m_mc_pow_emin;
		}
		else {
			m_mc_exponent   = 0.0;
			m_mc_pow_emin   = std::log(emin.MeV());
			m_mc_pow_ewidth = std::log(emax.MeV()) - m_mc_pow_emin;
		}

	}

    // Return
    return;
}
