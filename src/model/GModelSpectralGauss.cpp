/***************************************************************************
 *         GModelSpectralGauss.cpp - Spectral Gaussian model class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Christoph Deil & Ellis Owen                      *
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
 * @file GModelSpectralGauss.cpp
 * @brief Gaussian spectral model class implementation
 * @author Christoph Deil & Ellis Owen
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GModelSpectralGauss.hpp"
#include "GModelSpectralRegistry.hpp"
#include "GFunction.hpp"
#include "GIntegral.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpectralGauss    g_spectral_gauss_seed;
const GModelSpectralRegistry g_spectral_gauss_registry(&g_spectral_gauss_seed);

/* __ Method name definitions ____________________________________________ */
#define G_FLUX                "GModelSpectralGauss::flux(GEnergy&, GEnergy&)"
#define G_EFLUX              "GModelSpectralGauss::eflux(GEnergy&, GEnergy&)"
#define G_MC             "GModelSpectralGauss::mc(GEnergy&, GEnergy&, GRan&)"
#define G_READ                      "GModelSpectralGauss::read(GXmlElement&)"
#define G_WRITE                    "GModelSpectralGauss::write(GXmlElement&)"

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
GModelSpectralGauss::GModelSpectralGauss(void) : GModelSpectral()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs constant spectral model by extracting information from an XML
 * element. See the read() method for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelSpectralGauss::GModelSpectralGauss(const GXmlElement& xml) :
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
 * @brief Constructor
 *
 * @param[in] prefactor Power law pre factor (in ph/cm2/s/MeV).
 * @param[in] mean Mean energy.
 * @param[in] sigma Energy width.
 ***************************************************************************/
GModelSpectralGauss::GModelSpectralGauss(const double&  prefactor,
                                         const GEnergy& mean,
                                         const GEnergy& sigma) :
                     GModelSpectral()
{
    // Initialise members
    init_members();

    // Set parameters
    m_norm.value(prefactor);
    m_mean.value(mean.MeV());
    m_sigma.value(sigma.MeV());

    // Autoscale parameters
    autoscale();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Spectral Gaussian model.
 ***************************************************************************/
GModelSpectralGauss::GModelSpectralGauss(const GModelSpectralGauss& model) :
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
GModelSpectralGauss::~GModelSpectralGauss(void)
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
 * @param[in] model Gauss spectral model.
 * @return Gauss spectral model.
 ***************************************************************************/
GModelSpectralGauss& GModelSpectralGauss::operator=(const GModelSpectralGauss& model)
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
 * @brief Clear Gaussian spectral model
 ***************************************************************************/
void GModelSpectralGauss::clear(void)
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
 * @brief Clone Gaussian spectral model
 *
 * @return Pointer to deep copy of Gaussian spectral model.
 ***************************************************************************/
GModelSpectralGauss* GModelSpectralGauss::clone(void) const
{
    // Clone Gaussian spectral model
    return new GModelSpectralGauss(*this);
}


/***********************************************************************//**
 * @brief Evaluate model value
 *
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @return Model value (ph/cm2/s/MeV).
 *
 * Evaluates
 *
 * \f[
<<<<<<< HEAD
 * \frac{dN}{dE} = \frac{\norm}{\sqrt{2\pi}\sigma}
 *                 \exp(\frac{(E-\bar{E})^2}{2\sigma^2})
=======
 * \frac{dN}{dE}=\frac{\norm}{\sqrt{2\pi}\sigma}\exp(\frac{-(E-\bar{E})^2}{2\sigma^2})
>>>>>>> Implements eval_gradients in GModelSpectralGauss - not yet checked
 * \f]
 ***************************************************************************/
double GModelSpectralGauss::eval(const GEnergy& srcEng,
                                 const GTime&   srcTime) const
{
    double energy = srcEng.MeV();
    double norm   = m_norm.value();
    double mean   = m_mean.value();
    double sigma  = m_sigma.value();

    // Compute function value
    double term1 = (norm / sigma) * gammalib::inv_sqrt2pi;
    double term2 = (energy - mean) * (energy - mean) / (2.0 * sigma * sigma);
    double value = term1 * std::exp(-term2);

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate model value and gradient
 *
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @return Model value (ph/cm2/s/MeV).
 *
 * This method simply calls the eval() method as no analytical gradients will
 * be computed. See the eval() method for details.
 *
<<<<<<< HEAD
 * TODO: update docstring.
=======
 * TODO: Implement this
>>>>>>> efa5ea566f0bf64aca32f7e22a66761cbccf2f80
 ***************************************************************************/
double GModelSpectralGauss::eval_gradients(const GEnergy& srcEng,
                                           const GTime&   srcTime)
{
    	double energy = srcEng.MeV();
    	double norm = m_norm.value();
    	double mean = m_mean.value();
    	double sigma = m_sigma.value();

	// Update the evaluation cache
	    update_eval_cache(srcEng);

<<<<<<< HEAD
	    // Compute function terms
	    double term1 = (norm / sigma) * gammalib::inv_sqrt2pi;
	    double term2 = (1 / sigma) * gammalib::inv_sqrt2pi;
	    double term3 = (energy - mean) * (energy - mean) / (2 * sigma * sigma);
	    double term4 = (energy - mean) / (sigma * sigma);
	    double term5 = (norm / (sigma * sigma)) * gammalib::inv_sqrt2pi;

	    // Compute function value
	    double value = term1 * std::exp(- term3);


	    // Compute partial derivatives with respect to the parameter factor
	    // values (partial differentials were determined analytically).

	    double g_norm  = term2 * std::exp(- term3);
	    double g_mean  = term1 * term4 * std::exp(- term3);
	    double g_sigma = - term5 * std::exp(- term3) * (1 - (2 * term3));

		// Set gradients
=======
	    // Compute function value
	    double term1 = (norm / sigma) * gammalib::inv_sqrt2pi;
	    double term2 = (energy - mean) * (energy - mean) / (2 * sigma * sigma);
	    double value = term1 * std::exp(- term2);

	    // Compute partial derivatives with respect to the parameter factor
	    // values. The partial derivatives with respect to the parameter
	    // values are obtained by division by the scale factor.
	    double g_norm  = (m_norm.is_free())
	                     ? m_norm.scale() * m_last_power : 0.0;
	    double g_mean = (m_index.is_free())
	                     ? value * m_index.scale() * std::log(m_last_e_norm) : 0.0;
	    double g_sigma  = (m_ecut.is_free())
	                     ? value * m_last_e_cut / m_ecut.factor_value() : 0.0;

	    // Set gradients
>>>>>>> efa5ea566f0bf64aca32f7e22a66761cbccf2f80
	    m_norm.factor_gradient(g_norm);
	    m_mean.factor_gradient(g_mean);
	    m_sigma.factor_gradient(g_sigma);

	    // Compile option: Check for NaN/Inf
	    #if defined(G_NAN_CHECK)
	    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
<<<<<<< HEAD
	        std::cout << "*** ERROR: GModelSpectralGauss::eval_gradients";
=======
	        std::cout << "*** ERROR: GModelSpectralExpPlaw::eval_gradients";
>>>>>>> efa5ea566f0bf64aca32f7e22a66761cbccf2f80
	        std::cout << "(srcEng=" << srcEng;
	        std::cout << ", srcTime=" << srcTime << "):";
	        std::cout << " NaN/Inf encountered";
	        std::cout << " (value=" << value;
<<<<<<< HEAD
	        std::cout << ", norm=" << m_last_norm;
	        std::cout << ", mean=" << m_last_mean;
	        std::cout << ", Sigma=" << m_last_sigma;
=======
	        std::cout << ", e_norm=" << m_last_e_norm;
	        std::cout << ", e_cut=" << m_last_e_cut;
	        std::cout << ", power=" << m_last_power;
>>>>>>> efa5ea566f0bf64aca32f7e22a66761cbccf2f80
	        std::cout << ")" << std::endl;
	    }
	    #endif

	    // Return
	    return value;
	}

}


/***********************************************************************//**
 * @brief Returns model photon flux between [emin, emax] (ph/cm2/s)
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
double GModelSpectralGauss::flux(const GEnergy& emin,
                                 const GEnergy& emax) const
{
    // Initialise flux
    double flux = 0.0;
    
    // Compute only if integration range is valid
    if (emin < emax) {

        // Precomputations
        double energy_min = emin.MeV();
        double energy_max = emax.MeV();
        double norm       = m_norm.value();
        double mean       = m_mean.value();
        double sigma      = m_sigma.value();
        double denom      = 1.0 / (gammalib::sqrt_two*sigma);
        double zmin       = (energy_min - mean) * denom;
        double zmax       = (energy_max - mean) * denom;

        // Compute flux for a constant model
        flux = norm*(gammalib::erfcc(zmin) - gammalib::erfcc(zmax))/2.0;
    
    } // endif: integration range was valid

    // Return
    return flux;
}


/***********************************************************************//**
 * @brief Returns model energy flux between [emin, emax] (erg/cm2/s)
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
double GModelSpectralGauss::eflux(const GEnergy& emin,
                                  const GEnergy& emax) const
{
    // Initialise flux
    double eflux = 0.0;
    
    // Compute only if integration range is valid
    if (emin < emax) {

    	// Setup integration kernel
    	eflux_kernel integrand(m_norm.value(),  m_mean.value(),
    	                       m_sigma.value());
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
 * Returns Monte Carlo energy by randomly drawing from a constant between
 * the minimum and maximum photon energy.
 *
 * Method Used: Box-Muller transform, outlined here:
 * http://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
 *
 * Code from: http://www.design.caltech.edu/erik/Misc/Gaussian.html
 ***************************************************************************/
GEnergy GModelSpectralGauss::mc(const GEnergy& emin,
                                const GEnergy& emax,
                                const GTime&   time,
                                GRan&          ran) const
{
	double xmax = emax.MeV();
	double xmin = emin.MeV();
	double energy;

    // Throw an exception if energy range is invalid
    if (xmin >= xmax) {
        throw GException::erange_invalid(G_MC, xmin, xmax,
              "Minimum energy < maximum energy required.");
    }

    // Get uniform value between 0 and 1

    // Do ...
    do {
        double x1;
        double w;
    	do{
    		x1        = 2.0 * ran.uniform() - 1.0;
    		double x2 = 2.0 * ran.uniform() - 1.0;
    		w         = x1 * x1 + x2 * x2;
    	} while (w >= 1.0);

        // Do ...
        w = std::sqrt( (-2.0 * std::log( w ) ) / w );

        // Compute ...
    	double val = x1 * w;

    	// Map into [emin,emax] range
    	energy = m_sigma.value() * val + m_mean.value();

    } while ((xmin <= energy) && (energy < xmax));

    // Return energy
    return GEnergy(energy, "MeV");
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element containing Gaussian model information.
 *
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter name found in XML element.
 *
 * Read the spectral Gaussian information from an XML element.
 * The format of the XML elements is:
 *
 *     <spectrum type="Gaussian">
 *       <parameter name="Prefactor" scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Mean"      scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Sigma"     scale=".." value=".." min=".." max=".." free=".."/>
 *     </spectrum>
 *
 ***************************************************************************/
void GModelSpectralGauss::read(const GXmlElement& xml)
{
    // Set number of parameters
    const int n_pars = 3;

    // Verify that XML element has exactly 3 parameters
    if (xml.elements() != n_pars || xml.elements("parameter") != n_pars) {
        throw GException::model_invalid_parnum(G_READ, xml,
              "Gaussian model requires exactly 3 parameters.");
    }

    // Extract model parameters
    int npar[] = {0, 0, 0};
    for (int i = 0; i < n_pars; ++i) {

        // Get parameter element
        const GXmlElement* par = xml.element("parameter", i);

        // Handle prefactor
        if (par->attribute("name") == "Prefactor") {
            m_norm.read(*par);
            npar[0]++;
        }

        // Handle mean
        else if (par->attribute("name") == "Mean") {
            m_mean.read(*par);
            npar[1]++;
        }

        // Handle sigma
        else if (par->attribute("name") == "Sigma") {
            m_sigma.read(*par);
            npar[2]++;
        }

    } // endfor: looped over all parameters

    // Verify that all parameters were found
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1) {
        throw GException::model_invalid_parnames(G_READ, xml,
              "Require \"Prefactor\", \"Mean\", and \"Sigma\""
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
 *            Existing XML element is not of type "Gaussian"
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters or nodes found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Writes the spectral information into an XML element. The format of the XML
 * element is
 *
 *     <spectrum type="Gaussian">
 *       <parameter name="Prefactor" scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Mean"     scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Sigma"    scale=".." value=".." min=".." max=".." free=".."/>
 *     </spectrum>
 ***************************************************************************/
void GModelSpectralGauss::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "") {
        xml.attribute("type", "Gaussian");
    }

    // Verify model type
    if (xml.attribute("type") != "Gaussian") {
        throw GException::model_invalid_spectral(G_WRITE, xml.attribute("type"),
              "Spectral model is not of type \"Gaussian\".");
    }

    // If XML element has 0 nodes then append 3 parameter nodes
    if (xml.elements() == 0) {
        xml.append(GXmlElement("parameter name=\"Prefactor\""));
        xml.append(GXmlElement("parameter name=\"Mean\""));
        xml.append(GXmlElement("parameter name=\"Sigma\""));
    }

    // Verify that XML element has exactly 3 parameters
    if (xml.elements() != 3 || xml.elements("parameter") != 3) {
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Power law model requires exactly 3 parameters.");
    }

    // Set or update model parameter attributes
    int npar[] = {0, 0, 0};
    for (int i = 0; i < 3; ++i) {

        // Get parameter element
        GXmlElement* par = xml.element("parameter", i);

        // Handle prefactor
        if (par->attribute("name") == "Prefactor") {
            npar[0]++;
            m_norm.write(*par);
        }

        // Handle mean
        else if (par->attribute("name") == "Mean") {
            npar[1]++;
            m_mean.write(*par);
        }

        // Handle sigma
        else if (par->attribute("name") == "Sigma") {
            npar[2]++;
            m_sigma.write(*par);
        }


    } // endfor: looped over all parameters

    // Check of all required parameters are present
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1 ) {
        throw GException::model_invalid_parnames(G_WRITE, xml,
              "Require \"Prefactor\", \"Mean\" and \"Sigma\""
              " parameters.");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print spectral model information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing spectral model information.
 ***************************************************************************/
std::string GModelSpectralGauss::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpectralGauss ===");

        // Append model content
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
void GModelSpectralGauss::init_members(void)
{
    // Initialise normalisation
    m_norm.clear();
    m_norm.name("Prefactor");
    m_norm.unit("ph/cm2/s/MeV");
    m_norm.scale(1.0);
    m_norm.value(1.0);          // default: 1.0
    m_norm.min(0.0);            // min:     0.0
    m_norm.free();
    m_norm.gradient(0.0);
    m_norm.has_grad(false);

    // Initialise mean energy
    m_mean.clear();
    m_mean.name("Mean");
    m_mean.unit("MeV");
    m_mean.scale(1.0);
    m_mean.value(1000.0);       // default: 1000.0
    m_mean.min(0.1);            // min:     0.1
    m_mean.free();
    m_mean.gradient(0.0);
    m_mean.has_grad(false);

    // Initialise energy width
    m_sigma.clear();
    m_sigma.name("Sigma");
    m_sigma.unit("MeV");
    m_sigma.scale(1.0);
    m_sigma.value(1000.0);       // default: 1000.0
    m_sigma.min(0.1);            // min:     0.1
    m_sigma.free();
    m_sigma.gradient(0.0);
    m_sigma.has_grad(false);

    // Initialise eval cache
    m_last_norm = 0.0;
    m_last_mean  = 0.0;
    m_last_sigma  = 0.0;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);
    m_pars.push_back(&m_mean);
    m_pars.push_back(&m_sigma);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Spectral constant model.
 ***************************************************************************/
void GModelSpectralGauss::copy_members(const GModelSpectralGauss& model)
{
    // Copy members
    m_norm  = model.m_norm;
    m_mean  = model.m_mean;
    m_sigma = model.m_sigma;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);
    m_pars.push_back(&m_mean);
    m_pars.push_back(&m_sigma);

    // Copy eval cache
    m_last_norm = model.m_last_norm;
    m_last_mean = model.m_last_mean;
    m_last_sigma = model.m_last_sigma;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralGauss::free_members(void)
{
    // Return
    return;
}

/***********************************************************************//**
 * @brief Kernel for energy flux integration
 *
 * @param[in] energy Energy (MeV).
 ***************************************************************************/
double GModelSpectralGauss::eflux_kernel::eval(const double& energy)
{
    // Evaluate function value
    double term1 = (m_norm / m_sigma) * gammalib::inv_sqrt2pi;
    double term2 = (energy - m_mean) * (energy - m_mean) / (2.0 * m_sigma * m_sigma);
    double value = term1 * std::exp(- term2);

    // Return value
    return value;
}

/***********************************************************************//**
 * @brief Update eval precomputation cache
 *
 * @param[in] energy Energy.
 *
 * Updates the precomputation cache for eval() and eval_gradients() methods.
 *
<<<<<<< HEAD
=======
 * TODO integrate this to correct form
>>>>>>> efa5ea566f0bf64aca32f7e22a66761cbccf2f80
 ***************************************************************************/
void GModelSpectralGauss::update_eval_cache(const GEnergy& energy) const
{
    // Get parameter values (takes 3 multiplications which are difficult
    // to avoid)
<<<<<<< HEAD
    double norm = m_norm.value();
    double mean  = m_mean.value();
    double sigma = m_sigma.value();

    // If the energy or one of the parameters norm, mean or sigma
    // energy has changed then recompute the cache
    if ((m_last_norm != norm) ||
        (m_last_mean  != mean)  ||
        (m_last_sigma   != sigma)) {

        // Store actual energy and parameter values
        m_last_norm = norm;
        m_last_mean  = mean;
        m_last_sigma   = sigma;

        // Compute and store value
        double eng    = energy.MeV();
        m_last_norm = eng / m_last_norm;
        m_last_mean  = eng / m_last_mean;
        m_last_sigma  = eng/ m_last_sigma;
=======
    double index = m_index.value();
    double ecut  = m_ecut.value();
    double pivot = m_pivot.value();

    // If the energy or one of the parameters index, cut-off or pivot
    // energy has changed then recompute the cache
    if ((m_last_energy != energy) ||
        (m_last_index  != index)  ||
        (m_last_ecut   != ecut)   ||
        (m_last_pivot  != pivot)) {

        // Store actual energy and parameter values
        m_last_energy = energy;
        m_last_index  = index;
        m_last_ecut   = ecut;
        m_last_pivot  = pivot;

        // Compute and store value
        double eng    = energy.MeV();
        m_last_e_norm = eng / m_last_pivot;
        m_last_e_cut  = eng / m_last_ecut;
        m_last_power  = std::pow(m_last_e_norm, m_last_index) *
                        std::exp(-m_last_e_cut);
>>>>>>> efa5ea566f0bf64aca32f7e22a66761cbccf2f80

    } // endif: recomputation was required

    // Return
    return;
}


