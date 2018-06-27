/***************************************************************************
 *  GModelSpectralExponential.cpp - Exponential spectral model class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2018- by Lluigi Tibaldo                                  *
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
 * @file GModelSpectralExponential.cpp
 * @brief Exponential spectral model class implementation
 * @author Luigi Tibaldo
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
#include "GEnergies.hpp"
#include "GModelSpectralExponential.hpp"
#include "GModelSpectralNodes.hpp"
#include "GModelSpectralRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpectralExponential g_spectral_multi_seed;
const GModelSpectralRegistry    g_spectral_multi_registry(&g_spectral_multi_seed);

/* __ Method name definitions ____________________________________________ */
#define G_MC  "GModelSpectralExponential::mc(GEnergy&, GEnergy&, GTime&, "\
                                                                     "GRan&)"
#define G_WRITE           "GModelSpectralExponential::write(GXmlElement&)"
#define G_SETEXPONENT     "GModelSpectralExponential::exponent(GModelSpectral&)"
#define G_RETURNEXPONENT  "GModelSpectralExponential::exponent(void)"


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GModelSpectralExponential::GModelSpectralExponential(void) :
                              GModelSpectral()
{
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element containing spectral model information.
 *
 * Constructs a Exponential spectral model by extracting information from
 * an XML element. See the read() method for more information about the
 * expected structure of the XML element.
 ***************************************************************************/
GModelSpectralExponential::GModelSpectralExponential(const GXmlElement& xml) :
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
 * @param[in] model Exponential spectral model.
 ***************************************************************************/
GModelSpectralExponential::GModelSpectralExponential(const GModelSpectralExponential& model) :
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
GModelSpectralExponential::~GModelSpectralExponential(void)
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
 * @param[in] model Exponential spectral model.
 * @return Exponential spectral model.
 ***************************************************************************/
GModelSpectralExponential& GModelSpectralExponential::operator=(const GModelSpectralExponential& model)
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
 * @brief Clear Exponential spectral model
 ***************************************************************************/
void GModelSpectralExponential::clear(void)
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
 * @brief Clone Exponential spectral model model
 *
 * @return Pointer to deep copy of Exponential spectral model.
 ***************************************************************************/
GModelSpectralExponential* GModelSpectralExponential::clone(void) const
{
    // Clone spectral power law model
    return new GModelSpectralExponential(*this);
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
 *    \exp{M}(\rm srcEng, srcTime)
 * \f]
 *
 * where \f${M}\f$ is the exponent model component.
 *
 * If the @p gradients flag is true the method will also compute the partial
 * derivatives of each parameter with respect to the parameters using
 *
 * \f[
 *    \frac{\delta S}{\delta P_{\rm i}}\exp{M}
 * \f]
 *
 * where \f${P_{\rm i}}\f$ is the i-th parameter.
 ***************************************************************************/
double GModelSpectralExponential::eval(const GEnergy& srcEng,
                                          const GTime&   srcTime,
                                          const bool&    gradients) const
{
    // Initialise result
    double value = 0.0;

    // Check if exponent is defined
    if (m_exponent != NULL) {

        // Calculate exponent value
        value = m_exponent->eval(srcEng, srcTime, gradients);
        // calculate exponential
        value = exp(value);
        
		// Modify gradients if requested
		if (gradients) {
	
				// Loop over model parameters
				for (int ipar = 0; ipar < m_exponent->size(); ++ipar) {
	
					// Get reference to model parameter
					GModelPar& par = m_exponent->operator[](ipar);
	
					// Scale parameter gradient
					par.gradient(par.gradient()*value);
	
				} // endfor: loop over model parameters
				
		}// endif: compute grdients
	
	}//endif compute value and gradients

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GModelSpectralExponential::eval";
        std::cout << "(srcEng=" << srcEng;
        std::cout << ", srcTime=" << srcTime << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
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
 * Computes the photon flux of Exponential spectral model
 ***************************************************************************/
double GModelSpectralExponential::flux(const GEnergy& emin,
                                          const GEnergy& emax) const
{
    // Initialise flux
    double flux = 0.0;

    // Compute only if integration range is valid and exponent are available
    if (emin < emax && m_exponent != NULL) {

        // Initialise function to integrate
        flux_kern kernel(m_exponent);

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
 * Computes the energy flux of Exponential spectral model
 ***************************************************************************/
double GModelSpectralExponential::eflux(const GEnergy& emin,
                                           const GEnergy& emax) const
{
    // Initialise eflux
    double eflux = 0.0;

    // Compute only if integration range is valid and exponent are available
    if (emin < emax && m_exponent != NULL) {

        // Initialise function to integrate
        eflux_kern kernel(m_exponent);

        // Initialise integral class with function
        GIntegral integral(&kernel);

        // Set integration precision
        integral.eps(1.0e-8);

        // Calculate integral between emin and emax
        eflux = integral.romberg(emin.MeV(), emax.MeV());

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
 * Returns Monte Carlo energy by randomly drawing from a Exponential
 * spectral model.
 ***************************************************************************/
GEnergy GModelSpectralExponential::mc(const GEnergy& emin,
                                         const GEnergy& emax,
                                         const GTime&   time,
                                         GRan&          ran) const
{
    // Throw an exception if energy range is invalid
    if (emin >= emax) {
        throw GException::erange_invalid(G_MC, emin.MeV(), emax.MeV(),
              "Minimum energy < maximum energy required.");
    }

    // Throw exception if exponent is undefined
    if (m_exponent == 0) {
        std::string msg = "Exponent is undefined";
    	throw GException::runtime_error(G_MC, msg);
    }

    // Update MC cache
    update_mc_cache(emin, emax);

    // Set energy
    GEnergy energy = m_mc_spectrum.mc(emin, emax, time, ran);

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
void GModelSpectralExponential::read(const GXmlElement& xml)
{

	// Get exponent XML element
	const GXmlElement* spec = xml.element("spectrum",0);
	
    // Allocate a spectral registry object
    GModelSpectralRegistry registry;

    // Read spectral model
    GModelSpectral* ptr = registry.alloc(*spec);
    
    // Set spectral component as exponent
    exponent(*ptr);

    // Free spectral model
    delete ptr;
			
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
void GModelSpectralExponential::write(GXmlElement& xml) const
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
    
    // Check if exponent is defined
    if (m_exponent != NULL) {
    	
        // Create a spectrum node
        xml.append(GXmlElement("spectrum"));
        
        // Get spectrum node
        GXmlElement* spec = xml.element("spectrum", 0);
        
        // Create temporary copy of the exponent. This is a kludge to
        // write out the original parameters.
		GModelSpectral* cpy = m_exponent->clone();
		
		// Loop over parameters of model
		for (int i = 0; i < cpy->size(); ++i) {

			// Get model parameter and name
			GModelPar&  par     = (*cpy)[i];
			std::string parname = par.name();

			// Check if name contains colon
			if (gammalib::contains(parname, ":")) {

				// Split at the colon
				std::vector<std::string> splits = gammalib::split(parname, ":");

				// Use second part of the string to recover original
				// parameter name
				par.name(splits[1]);
			}

		} // endfor: loop over parameters

		// Write spectral component
		cpy->write(*spec);
		
        // Remove temporary copy
        delete cpy;
    	
    }

    // Return
    return;
}

/***********************************************************************//**
 * @brief Set exponent
 *
 * @param[in] spec Spectral model to use as exponent.
 *
 * Set a spectral model as exponent
 ***************************************************************************/
void GModelSpectralExponential::exponent(const GModelSpectral& spec)
{
	// Set exponent
	m_exponent = spec.clone();
	
	// Get number of spectral parameters from model
	int npars = m_exponent->size();
	
    // Store pointers to spectral parameters
    m_pars.clear();
	for (int ipar = 0; ipar < spec->size(); ++ipar) {

		// Get model parameter reference
		GModelPar& par = spec->operator[](ipar);

		// Append model parameter pointer to internal container
		m_pars.push_back(&par);
	}

    // Return
    return;
	
}


/***********************************************************************//**
 * @brief Return exponent
 * 
 * @return Pointer to spectral model.
 *
 * Returns the exponent of the Exponential spectral model.
 ***************************************************************************/
const GModelSpectral* GModelSpectralExponential::component(void) const
{
    // Return exponent
    return m_exponent;
}


/***********************************************************************//**
 * @brief Print Exponential spectral model information
 *
 * @param[in] chatter Chattiness.
 * @return String containing model information.
 ***************************************************************************/
std::string GModelSpectralExponential::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpectralExponential ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of parameters"));
        result.append(gammalib::str(size()));

        // Print parameter information
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
void GModelSpectralExponential::init_members(void)
{
    // Initialise model type
    m_type = "Exponential";

    // Clear exponent
    m_exponent.clear();

    // Clear MC cache
    m_mc_spectrum.clear();
    m_mc_emin.clear();
    m_mc_emax.clear();
    m_mc_values.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Exponential spectral model.
 ***************************************************************************/
void GModelSpectralExponential::copy_members(const GModelSpectralExponential& model)
{
    // Copy members
    m_type        = model.m_type;

    // Copy MC cache
    m_mc_spectrum = model.m_mc_spectrum;
    m_mc_emin     = model.m_mc_emin;
    m_mc_emax     = model.m_mc_emax;
    m_mc_values   = model.m_mc_values;

    // Set exponent
    m_exponent.clear();
    exponent(model.m_exponent);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralExponential::free_members(void)
{
	// Free exponent
	if (m_exponent != NULL) delete m_exponent;
	
	// Signal free pointer
	m_exponent = NULL;
	
    // Return
    return;
	
}

/***********************************************************************//**
 * @brief Update Monte Carlo pre computation cache
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 *
 * Updates the precomputation cache for Monte Carlo simulations. In case that
 * the energy boundaries have changed or at least one of the model parameters
 * has changed the method computes a spectral node function which has 100
 * nodes per decade containing the Exponential spectral model values and
 * stores that into a Monte Carlo cache. This cache is then used by the mc()
 * method for simulations.
 ***************************************************************************/
void GModelSpectralExponential::update_mc_cache(const GEnergy& emin,
                                                   const GEnergy& emax) const

{
    // Check if one of the parameters has changed. If the dimension of the
    // parameter value cache differs from the number of parameters then notify
    // a change. This will clear the cache and store the parameter values
    // later
    bool par_changed = (m_mc_values.size() != size());
    if (par_changed == false) {
        for (int i = 0; i < size(); ++i) {
            if (m_pars[i]->value() != m_mc_values[i]) {
                par_changed = true;
                break;
            }
        }
    }

    // Update cache if energy range or parameters have changed
    if (par_changed || emin != m_mc_emin || emax != m_mc_emax) {

        // Store energy range
        m_mc_emin = emin;
        m_mc_emax = emax;

        // If parameters have changed then store the current parameter
        // values for a comparison check for the next method call
        if (par_changed) {
            m_mc_values.clear();
            for (int i = 0; i < size(); ++i) {
                m_mc_values.push_back(m_pars[i]->value());
            }
        }

        // Clear spectral nodes
        m_mc_spectrum.clear();

        // Compute number of nodes. We use here 100 nodes per log energy and
        // assure that at least 100 nodes are used.
        int nodes = int((emax.log10MeV() - emin.log10MeV()) * 100.0);
        if (nodes < 100) {
            nodes = 100;
        }

        // Initialise energy array with fixed number of nodes
        GEnergies energies = GEnergies(nodes, m_mc_emin, m_mc_emax, true);

        // Append nodes to spectral function
        for (int i = 0; i < energies.size(); ++i) {
            m_mc_spectrum.append(energies[i], eval(energies[i]));
        }

    } // endif: emin and emax have changed

    // Return
    return;
}
