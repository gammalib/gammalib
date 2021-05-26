/***************************************************************************
 *         GModelSpectralComposite.cpp - Spectral composite model class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016-2021 by Michael Mayer                               *
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
 * @file GModelSpectralComposite.cpp
 * @brief Composite spectral model class implementation
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
#include "GModelSpectralComposite.hpp"
#include "GModelSpectralRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpectralComposite g_spectral_comp_seed;
const GModelSpectralRegistry  g_spectral_comp_registry(&g_spectral_comp_seed);

/* __ Method name definitions ____________________________________________ */
#define G_MC "GModelSpectralComposite::mc(GEnergy&, GEnergy&, GTime&, GRan&)"
#define G_WRITE                "GModelSpectralComposite::write(GXmlElement&)"
#define G_COMPONENT_INDEX          "GModelSpectralComposite::component(int&)"
#define G_COMPONENT_NAME   "GModelSpectralComposite::component(std::string&)"
#define G_APPEND          "GModelSpectralComposite::append(GModelSpectral&, "\
                                                              "std::string&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_UPDATE_MC_CACHE                  //!< Debug MC cache update


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GModelSpectralComposite::GModelSpectralComposite(void) : GModelSpectral()
{
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element containing model information.
 *
 * Constructs a composite spectral model by extracting information from an
 * XML element.
 ***************************************************************************/
GModelSpectralComposite::GModelSpectralComposite(const GXmlElement& xml) :
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
 * @param[in] model Composite spectral model.
 ***************************************************************************/
GModelSpectralComposite::GModelSpectralComposite(const GModelSpectralComposite& model) :
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
GModelSpectralComposite::~GModelSpectralComposite(void)
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
 * @param[in] model Composite spectral model.
 * @return Composite spectral model.
 ***************************************************************************/
GModelSpectralComposite& GModelSpectralComposite::operator=(const GModelSpectralComposite& model)
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
 * @brief Clear composite spectral model
 ***************************************************************************/
void GModelSpectralComposite::clear(void)
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
 * @brief Clone composite spectral model
 *
 * @return Pointer to deep copy of composite spectral model.
 ***************************************************************************/
GModelSpectralComposite* GModelSpectralComposite::clone(void) const
{
    // Clone spectral power law model
    return new GModelSpectralComposite(*this);
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
 *    \sum_{i=0}^{N} {M_{\rm i}}(\rm srcEng, srcTime)
 * \f]
 *
 * where \f${M_{\rm i}}\f$ is the model evaluation of the i-th model
 *
 * If the @p gradients flag is true the method will also compute the
 * partial derivatives of each model component with respect to the
 * parameters.
 ***************************************************************************/
double GModelSpectralComposite::eval(const GEnergy& srcEng,
                                     const GTime&   srcTime,
                                     const bool&    gradients) const
{
    // Initialise result
	double value = 0.0;

	// Loop over model components
	for (int i = 0; i < m_spectral.size(); ++i) {

		// Add model component
		value += m_spectral[i]->eval(srcEng, srcTime, gradients);

	} // endfor: loop over model components

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GModelSpectralComposite::eval";
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
 * Computes the sum of the photon fluxes in all components.
 ***************************************************************************/
double GModelSpectralComposite::flux(const GEnergy& emin,
                                     const GEnergy& emax) const
{
    // Initialise flux
    double flux = 0.0;

    // Compute only if integration range is valid
    if (emin < emax) {

    	// Loop over model components
    	for (int i = 0; i < m_spectral.size(); ++i) {

    		// Add model component
    		flux += m_spectral[i]->flux(emin, emax);

    	} // endfor: loop over model components

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
 * Computes the sum of the energy fluxes in all components.
 ***************************************************************************/
double GModelSpectralComposite::eflux(const GEnergy& emin,
                                      const GEnergy& emax) const
{
   // Initialise eflux
	double eflux = 0.0;

	// Compute only if integration range is valid
	if (emin < emax) {

		// Loop over model components
		for (int i = 0; i < m_spectral.size(); ++i) {

			// Add model component
			eflux += m_spectral[i]->flux(emin, emax);

		} // endfor: loop over model components

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
 * Returns Monte Carlo energy by randomly drawing from a composite spectral
 * model.
 ***************************************************************************/
GEnergy GModelSpectralComposite::mc(const GEnergy& emin,
                                    const GEnergy& emax,
                                    const GTime&   time,
                                    GRan&          ran) const
{
    // Throw an exception if energy range is invalid
    if (emin >= emax) {
        throw GException::erange_invalid(G_MC, emin.MeV(), emax.MeV(),
              "Minimum energy < maximum energy required.");
    }

    // Throw exception if model container is empty
    if (m_spectral.size() == 0) {
    	throw GException::runtime_error(G_MC,
    	              "Composite spectral model is empty. It is required to have at"\
					  " least one spectral model to run simulations");
    }

    // Update MC cache
    update_mc_cache(emin, emax);

    // Get random value
    double u = ran.uniform();

    // Initialise index of model to use for simulation
    int index = 0;

    // Loop over spectral components and compute relative probabilites
    for (int i = 0; i < m_mc_probs.size(); ++i) {
    	if (u <= m_mc_probs[i]) {
    		index = i;
    		break;
    	}
    }

    // Set energy
    GEnergy energy = m_spectral[index]->mc(emin, emax, time, ran);

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
void GModelSpectralComposite::read(const GXmlElement& xml)
{
    // Get number of spectral components
    int n_spectrals = xml.elements("spectrum");

    // Loop over spectral elements
    for (int i = 0; i < n_spectrals; ++i) {

        // Get spectral XML element
        const GXmlElement* spec = xml.element("spectrum", i);

        // Initialise a spectral registry object
        GModelSpectralRegistry registry;

        // Read spectral model
        GModelSpectral* ptr = registry.alloc(*spec);

        // Get component attribute from XML file
        std::string component_name = spec->attribute("component");

        // Append spectral component to container
        append(*ptr, component_name);

        // Free spatial model
        delete ptr;

    } // endfor: loop over components

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
void GModelSpectralComposite::write(GXmlElement& xml) const
{
    // Verify model type
    gammalib::xml_check_type(G_WRITE, xml, type());

    // Loop over model components
    for (int i = 0; i < m_spectral.size(); i++) {

        // Write spectral model
        if (m_spectral[i] != NULL) {

            // Create new spectrum node
            xml.append(GXmlElement("spectrum"));

            // Get new spectrum node
            GXmlElement* spec = xml.element("spectrum", xml.elements("spectrum")-1);

            // Create temporary copy of the spectral model
            // This is a kluge to write out the original parameters.
            GModelSpectral* cpy = m_spectral[i]->clone();

            // Loop over parameters of model
            for (int j = 0; j < cpy->size(); ++j) {

                // Get model parameter and name
                GModelPar& par = (*cpy)[j];
                std::string parname = par.name();

                // Check if name contains colon
                if (gammalib::contains(parname, ":")) {

                    // Split at the colon
                    std::vector<std::string> splits = gammalib::split(parname, ":");

                    // Use second part of the string to recover original parameter name
                    par.name(splits[1]);
                }

            } // endfor: loop over parameters

            // Write spectral component
            cpy->write(*spec);

            // Add component name if previously available
            if (m_components[i] != gammalib::str(i+1)) {
                spec->attribute("component", m_components[i]);
            }

            // Remove temporary copy
            delete cpy;

        } // endif: spectral model was not NULL

    } // endfor: loop over model components

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append spectral component
 *
 * @param[in] spec Spectral model component.
 * @param[in] name Name of spectral component.
 *
 * Appends a spectral component to the composite model
 ***************************************************************************/
void GModelSpectralComposite::append(const GModelSpectral& spec,
                                     const std::string&    name)
{
    // Append model container
    m_spectral.push_back(spec.clone());

    // Get index of latest model
    int index = m_spectral.size()-1;

    // Use model index if component name is empty
    std::string component_name = !name.empty() ? name
                                               : gammalib::str(m_spectral.size());

    // Check if component name is unique, throw exception if not
    if (gammalib::contains(m_components, component_name)) {
    	std::string msg = "Attempt to append component \""+component_name+"\" "
                          "to composite spectral model, but a component with "
                          "the same name exists already. Each component needs "
                          "a unique name.";
		throw GException::invalid_value(G_APPEND, msg);
    }

    // Add component name (for now simple number)
    m_components.push_back(component_name);

    // Get number of spectral parameters from model
    int npars = m_spectral[index]->size();

    // Loop over model parameters
    for (int ipar = 0; ipar < npars; ++ipar) {

        // Get model parameter
        GModelPar* par = &(m_spectral[index]->operator[](ipar));

        // Modify parameter name
        par->name(component_name+":"+par->name());

        // Append model parameter with new name to internal container
        m_pars.push_back(par);

    } // endfor: loop over model parameters

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns spectral component element
 *
 * @param[in] index Index of spectral component.
 * @return Spectral model.
 *
 * @exception GException::out_of_range
 *            Index is out of range.
 *
 * Returns a spectral component to the composite model
 ***************************************************************************/
const GModelSpectral* GModelSpectralComposite::component(const int& index) const
{
    // Check if index is in validity range
    if (index >= m_spectral.size() || index < 0) {
        throw GException::out_of_range(G_COMPONENT_INDEX, "Component Index",
                                       index, m_spectral.size());
    }

    // Return spectral component
    return m_spectral[index];
}


/***********************************************************************//**
 * @brief Returns pointer to specific spectral component
 *
 * @param[in] name Name of spectral component.
 * @return Spectral model.
 *
 * @exception GException::invalid_argument
 *            Spectral model not found.
 *
 * Returns a spectral component of the composite model
 ***************************************************************************/
const GModelSpectral* GModelSpectralComposite::component(const std::string& name) const
{
    // Check if model name is found
    int index = -1;
    for (int i = 0; i < m_components.size(); ++i) {
        if (m_components[i] == name) {
            index = i;
            break;
        }
    }

    // Check if component name was found
    if (index == -1) {
        std::string msg = "Model component \""+name+"\" not found in composite "
                          "spectral model.";
        throw GException::invalid_argument(G_COMPONENT_NAME, msg);
    }

    // Return spectral component
    return m_spectral[index];

}


/***********************************************************************//**
 * @brief Print composite spectral model information
 *
 * @param[in] chatter Chattiness.
 * @return String containing composite spectral model information.
 ***************************************************************************/
std::string GModelSpectralComposite::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

    	// Append header
        result.append("=== GModelSpectralComposite ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of components"));
        result.append(gammalib::str(components()));
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
void GModelSpectralComposite::init_members(void)
{
    // Initialise model type
    m_type = "Composite";

    // Clear spectral models
    m_spectral.clear();
    m_components.clear();

    // Clear MC cache
    m_mc_probs.clear();
    m_mc_emin.clear();
    m_mc_emax.clear();
    m_mc_values.clear();
    m_mc_flux = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model GModelSpectralComposite members which should be copied.
 ***************************************************************************/
void GModelSpectralComposite::copy_members(const GModelSpectralComposite& model)
{
    // Copy members
    m_type       = model.m_type;
    m_components = model.m_components;

    // Copy MC cache
    m_mc_probs   = model.m_mc_probs;
    m_mc_emin    = model.m_mc_emin;
    m_mc_emax    = model.m_mc_emax;
    m_mc_values  = model.m_mc_values;
    m_mc_flux    = model.m_mc_flux;

    // Copy pointer(s) of spectral component
    m_spectral.clear();

    // Clear parameters and copy the pointers from the clone
    m_pars.clear();
    for (int i = 0; i < model.components(); ++i) {

        // Clone spectral component
        m_spectral.push_back(model.m_spectral[i]->clone());

    	// Retrieve spectral model
    	GModelSpectral* spec = m_spectral[i];

    	// Loop over parameters and store pointers
		for (int ipar = 0; ipar < spec->size(); ++ipar) {

            // Get model parameter reference
            GModelPar& par = spec->operator[](ipar);

            // Append model parameter pointer to internal container
            m_pars.push_back(&par);

        } // endfor: loop over parameters

    } // endfor: loop over spectral models

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralComposite::free_members(void)
{
    // Free memory
    for (int i = 0; i < m_spectral.size(); ++i) {

        // Delete component i
        if (m_spectral[i] != NULL) {
            delete m_spectral[i];
        }

        // Signal free pointer
        m_spectral[i] = NULL;
	}

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
void GModelSpectralComposite::update_mc_cache(const GEnergy& emin,
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

        // Compute MC flux
        m_mc_flux = flux(m_mc_emin, m_mc_emax);

        // Initialise sum and probabilites
	    double sum = 0.0;
	    m_mc_probs.clear();

	    // Loop over spectral components and compute relative probabilites
	    for (int i = 0; i < m_spectral.size(); ++i) {

	    	// Relative probability
	    	double prob = m_spectral[i]->flux(emin, emax) / m_mc_flux;

	    	// Add probability
	    	m_mc_probs.push_back(prob + sum);

	    	// Increment sum
	    	sum += prob;

	    } //endfor: looped over spectral components

        // Debug option: signal update
        #if defined(G_DEBUG_UPDATE_MC_CACHE)
        std::cout << "GModelSpectralComposite::update_mc_cache(";
        std::cout << emin.print() << "," << emax.print() << "):";
        std::cout << " flux=" << m_mc_flux;
        for (int i = 0; i < m_spectral.size(); ++i) {
            std::cout << " prob[" << i << "]=" << m_mc_probs[i];
        }
        std::cout << std::endl;
        #endif

	} // endif: emin and emax have changed

    // Return
    return;
}
