/***************************************************************************
 *         GModelSpectralComposite.cpp - Spectral composite model class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Michael Mayer                                    *
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
const GModelSpectralComposite     g_spectral_comp_seed;
const GModelSpectralRegistry g_spectral_comp_registry(&g_spectral_comp_seed);

/* __ Method name definitions ____________________________________________ */
#define G_MC "GModelSpectralComposite::mc(GEnergy&, GEnergy&, GTime&, GRan&)"
#define G_READ                  "GModelSpectralComposite::read(GXmlElement&)"
#define G_WRITE                "GModelSpectralComposite::write(GXmlElement&)"

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
 * @param[in] xml XML element containing position information.
 *
 * Constructs a power law spectral model by extracting information from an
 * XML element. See the read() method for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelSpectralComposite::GModelSpectralComposite(const GXmlElement& xml) : GModelSpectral()
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
GModelSpectralComposite::GModelSpectralComposite(const GModelSpectralComposite& model)
                                       : GModelSpectral(model)
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
 * @param[in] model Spectral power law model.
 * @return Spectral power law model.
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
 * @brief Clear spectral power law model
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
 * @brief Clone spectral power law model
 *
 * @return Pointer to deep copy of spectral power law model.
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
 * partial derivatives of each model component with respect to the parameters using
 *
 * @todo The method expects that energy!=0. Otherwise Inf or NaN may result.
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
 * Computes the sum of all photon fluxes of each individual component
 *
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
 * Computes the sum of all energy fluxes of each individual component
 *
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
 * Returns Monte Carlo energy by randomly drawing from a power law.
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

    // Set energy
    GEnergy energy;

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
	for(int i = 0; i < n_spectrals; ++i) {

		// Get spectral XML element
		const GXmlElement* spec = xml.element("spectrum", i);

		// Add spectral component
		add_component(*spec);

	} // endfor: loop over components

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
void GModelSpectralComposite::write(GXmlElement& xml) const
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

    // Loop over model components
    for(int i = 0; i < m_spectral.size(); i++) {

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
			std::cout<<m_components[i]<<" " <<gammalib::str(i+1)<<std::endl;
			if (m_components[i] != gammalib::str(i+1)) {
				spec->attribute("component", m_components[i]);
			}

			// Remove temporary copy
			delete cpy;
		}

    } // endfor: loop over model components

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print powerlaw information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing model information.
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

    // Copy pointer(s) of spectral component
    m_spectral.clear();
    for(int i = 0; i < model.components(); ++i) {
    	GModelSpectral* spec = (model.m_spectral[i] != NULL) ? model.m_spectral[i]->clone() : NULL;
    	m_spectral.push_back(spec);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralComposite::free_members(void)
{
    // Free memory
	for(int i = 0; i < m_spectral.size(); ++i) {

		// Delete component i
		if (m_spectral[i] != NULL) delete m_spectral[i];

		// Signal free pointer
		m_spectral[i] = NULL;
	}

    // Return
    return;
}


/***********************************************************************//**
 * @brief Add spectral component from XML element
 ***************************************************************************/
void GModelSpectralComposite::add_component(const GXmlElement& spec)
{
	// Initialise a spectral registry object
    GModelSpectralRegistry registry;

	// Read spectral model
	GModelSpectral* ptr = registry.alloc(spec);

	// Get component attribute from XML file
	std::string component_name = spec.attribute("component");

	// Append spectral component to container
	append(*ptr, component_name);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append spectral component
 *
 * Appends a spectral component to the composite model
 ***************************************************************************/
void GModelSpectralComposite::append(const GModelSpectral& spec, const std::string& name)
{
    // Append model container
    m_spectral.push_back(spec.clone());

    // Get index of latest model
    int index = m_spectral.size()-1;

    // Use model index if component name is empty
    std::string component_name = !name.empty() ? name : gammalib::str(m_spectral.size());

	// Add component name (for now simple number)
	m_components.push_back(component_name);

	// Get number of spectral parameters from model
	int npars = m_spectral[index]->size();

	// Loop over model parameters
	for (int i = 0; i < npars; ++i) {

		// Get model parameter
		GModelPar& par = (*m_spectral[index])[i];

		// Modify parameter name
		par.name(name+":"+par.name());

		// Append model parameter with new name to internal container
		m_pars.push_back(&par);

	} // endfor: loop over model parameters

}
