/***************************************************************************
 *  GModelSpectralMultiplicative.cpp - Multiplicative spectral model class *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016-2020 by Michael Mayer                               *
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
 * @file GModelSpectralMultiplicative.cpp
 * @brief Multiplicative spectral model class implementation
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
#include "GEnergies.hpp"
#include "GModelSpectralMultiplicative.hpp"
#include "GModelSpectralNodes.hpp"
#include "GModelSpectralRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpectralMultiplicative g_spectral_multi_seed;
const GModelSpectralRegistry       g_spectral_multi_registry(&g_spectral_multi_seed);

/* __ Method name definitions ____________________________________________ */
#define G_MC  "GModelSpectralMultiplicative::mc(GEnergy&, GEnergy&, GTime&, "\
                                                                     "GRan&)"
#define G_WRITE           "GModelSpectralMultiplicative::write(GXmlElement&)"
#define G_COMPONENT_INDEX     "GModelSpectralMultiplicative::component(int&)"
#define G_COMPONENT_NAME            "GModelSpectralMultiplicative::component"\
                                                             "(std::string&)"
#define G_APPEND     "GModelSpectralMultiplicative::append(GModelSpectral&, "\
                                                              "std::string&)"

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
GModelSpectralMultiplicative::GModelSpectralMultiplicative(void) :
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
 * Constructs a multiplicative spectral model by extracting information from
 * an XML element. See the read() method for more information about the
 * expected structure of the XML element.
 ***************************************************************************/
GModelSpectralMultiplicative::GModelSpectralMultiplicative(const GXmlElement& xml) :
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
 * @param[in] model Multiplicative spectral model.
 ***************************************************************************/
GModelSpectralMultiplicative::GModelSpectralMultiplicative(const GModelSpectralMultiplicative& model) :
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
GModelSpectralMultiplicative::~GModelSpectralMultiplicative(void)
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
 * @param[in] model Multiplicative spectral model.
 * @return Multiplicative spectral model.
 ***************************************************************************/
GModelSpectralMultiplicative& GModelSpectralMultiplicative::operator=(const GModelSpectralMultiplicative& model)
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
 * @brief Clear multiplicative spectral model
 ***************************************************************************/
void GModelSpectralMultiplicative::clear(void)
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
 * @brief Clone multiplicative spectral model model
 *
 * @return Pointer to deep copy of multiplicative spectral model.
 ***************************************************************************/
GModelSpectralMultiplicative* GModelSpectralMultiplicative::clone(void) const
{
    // Clone spectral power law model
    return new GModelSpectralMultiplicative(*this);
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
 *    \prod_{i=0}^{N} {M_{\rm i}}(\rm srcEng, srcTime)
 * \f]
 *
 * where \f${M_{\rm i}}\f$ is the i-th model component.
 *
 * If the @p gradients flag is true the method will also compute the partial
 * derivatives of each parameter of eachmodel component with respect to the
 * parameters using
 *
 * \f[
 *    \frac{\delta S}{\delta P_{\rm ij}}\prod_{\rm k\neq \rm i}^{n} M_{\rm k}
 * \f]
 *
 * where \f${P_{\rm ij}}\f$ is the j-th parameter of the i-th multiplicative
 * component, while \f${M_{\rm k}}\f$ is the k-th model component and n the
 * number of model components.
 *
 * @todo The method expects that energy!=0. Otherwise Inf or NaN may result.
 ***************************************************************************/
double GModelSpectralMultiplicative::eval(const GEnergy& srcEng,
                                          const GTime&   srcTime,
                                          const bool&    gradients) const
{
    // Initialise result
    double value = 0.0;

    // Check for available spectral components
    if (m_spectral.size()) {

        // Set first spectral component
        value = m_spectral[0]->eval(srcEng, srcTime, gradients);

        // Loop over model components
        for (int i = 1; i < m_spectral.size(); ++i) {
            value *= m_spectral[i]->eval(srcEng, srcTime, gradients);
        }

    } // endfor: loop over model components

    // Modify gradients if requested
    if (gradients) {

        // Loop over model components
        for (int i = 0; i < m_spectral.size(); ++i) {

            // Initialise scaling factor
            double factor = 1.0;

            // Loop over other model components and compute factor
            for (int j = 0; j < m_spectral.size(); ++j) {
                if (i != j) {
                    factor *= m_spectral[j]->eval(srcEng, srcTime, false);
                }
            }

            // Loop over model parameters
            for (int ipar = 0; ipar < m_spectral[i]->size(); ++ipar) {

                // Get reference to model parameter
                GModelPar& par = m_spectral[i]->operator[](ipar);

                // Scale parameter gradient
                par.gradient(par.gradient()*factor);

            } // endfor: loop over model parameters

        } // endfor: loop over models

    } //endif: compute gradients

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GModelSpectralMultiplicative::eval";
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
 * Computes the photon flux of multiplicative spectral model
 ***************************************************************************/
double GModelSpectralMultiplicative::flux(const GEnergy& emin,
                                          const GEnergy& emax) const
{
    // Initialise flux
    double flux = 0.0;

    // Compute only if integration range is valid and components are available
    if (emin < emax && m_spectral.size()) {

        // Initialise function to integrate
        flux_kern kernel(m_spectral);

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
 * Computes the energy flux of multiplicative spectral model
 ***************************************************************************/
double GModelSpectralMultiplicative::eflux(const GEnergy& emin,
                                           const GEnergy& emax) const
{
    // Initialise eflux
    double eflux = 0.0;

    // Compute only if integration range is valid and components are available
    if (emin < emax && m_spectral.size()) {

        // Initialise function to integrate
        eflux_kern kernel(m_spectral);

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
 * Returns Monte Carlo energy by randomly drawing from a multiplicative
 * spectral model.
 ***************************************************************************/
GEnergy GModelSpectralMultiplicative::mc(const GEnergy& emin,
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
        std::string msg = "Multiplicative spectral model is empty. At least"
                          " one spectral model is required for simulations.";
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
void GModelSpectralMultiplicative::read(const GXmlElement& xml)
{
    // Get number of spectral components
    int n_spectrals = xml.elements("spectrum");

    // Loop over spectral elements
    for (int i = 0; i < n_spectrals; ++i) {

        // Get spectral XML element
        const GXmlElement* spec = xml.element("spectrum", i);

        // Allocate a spectral registry object
        GModelSpectralRegistry registry;

        // Read spectral model
        GModelSpectral* ptr = registry.alloc(*spec);

        // Get component attribute from XML file
        std::string component_name = spec->attribute("component");

        // Append spectral component to container
        append(*ptr, component_name);

        // Free spectral model
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
 * @exception GException::model_invalid_spectral
 *            Existing XML element is not of the expected type.
 *
 * Writes the spectral information into an XML element.
 ***************************************************************************/
void GModelSpectralMultiplicative::write(GXmlElement& xml) const
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
    for (int i = 0; i < m_spectral.size(); i++) {

        // Write spectral model
        if (m_spectral[i] != NULL) {

            // Create new spectrum node
            xml.append(GXmlElement("spectrum"));

            // Get new spectrum node
            GXmlElement* spec = xml.element("spectrum", xml.elements("spectrum")-1);

            // Create temporary copy of the spectral model. This is a kluge to
            // write out the original parameters.
            GModelSpectral* cpy = m_spectral[i]->clone();

            // Loop over parameters of model
            for (int j = 0; j < cpy->size(); ++j) {

                // Get model parameter and name
                GModelPar&  par     = (*cpy)[j];
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
 * @param[in] name Name of spectral component (can be empty).
 *
 * @exception GException::invalid_value
 *            Invalid component name specified
 *
 * Appends a spectral component to the Multiplicative model
 ***************************************************************************/
void GModelSpectralMultiplicative::append(const GModelSpectral& spec,
                                          const std::string&    name)
{
    // Append model container
    m_spectral.push_back(spec.clone());

    // Get index of latest model
    int index = m_spectral.size()-1;

    // Use model index as component name if component name is empty
    std::string component_name = !name.empty() ? name
                                               : gammalib::str(m_spectral.size());

    // Check if component name is unique, throw exception if not
    if (gammalib::contains(m_components, component_name)) {
        std::string msg = "Attempt to append component with name \""+
                          component_name+"\" to multiplicative spectral model "
                          "container, but a component with the same name exists "
                          "already. Every component in the container needs a "
                          "unique name. On default the system will increment "
                          "an integer if no component name is provided.";
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

        // Prepend component name to parameter name
        par->name(component_name+":"+par->name());

        // Append model parameter with new name to internal container
        m_pars.push_back(par);

    } // endfor: loop over model parameters

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return spectral model component by index
 *
 * @param[in] index Index of spectral component.
 * @return Pointer to spectral model.
 *
 * Returns a component of the multiplicative spectral model by @p index.
 ***************************************************************************/
const GModelSpectral* GModelSpectralMultiplicative::component(const int& index) const
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
 * @brief Return spectral model component by name
 *
 * @param[in] name Name of spectral component.
 * @return Pointer to spectral model.
 *
 * Returns a component of the multiplicative spectral model by @p name.
 ***************************************************************************/
const GModelSpectral* GModelSpectralMultiplicative::component(const std::string& name) const
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
        throw GException::model_not_found(G_COMPONENT_NAME, name);
    }

    // Return spectral component
    return m_spectral[index];
}


/***********************************************************************//**
 * @brief Print multiplicative spectral model information
 *
 * @param[in] chatter Chattiness.
 * @return String containing model information.
 ***************************************************************************/
std::string GModelSpectralMultiplicative::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpectralMultiplicative ===");

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
void GModelSpectralMultiplicative::init_members(void)
{
    // Initialise model type
    m_type = "Multiplicative";

    // Clear spectral models
    m_spectral.clear();
    m_components.clear();

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
 * @param[in] model Multiplicative spectral model.
 ***************************************************************************/
void GModelSpectralMultiplicative::copy_members(const GModelSpectralMultiplicative& model)
{
    // Copy members
    m_type        = model.m_type;
    m_components  = model.m_components;

    // Copy MC cache
    m_mc_spectrum = model.m_mc_spectrum;
    m_mc_emin     = model.m_mc_emin;
    m_mc_emax     = model.m_mc_emax;
    m_mc_values   = model.m_mc_values;

    // Copy pointer(s) of spectral component
    m_spectral.clear();
    for (int i = 0; i < model.components(); ++i) {
    	m_spectral.push_back(model.m_spectral[i]->clone());
    }

    // Store pointers to spectral parameters
    m_pars.clear();
    for (int i = 0; i < model.components(); ++i) {

        // Retrieve spectral model
        GModelSpectral* spec = m_spectral[i];

        // Loop over parameters
        for (int ipar = 0; ipar < spec->size(); ++ipar) {

            // Get model parameter reference
            GModelPar& par = spec->operator[](ipar);

            // Append model parameter pointer to internal container
            m_pars.push_back(&par);
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralMultiplicative::free_members(void)
{
    // Free memory
    for (int i = 0; i < m_spectral.size(); ++i) {

        // Delete component i
        if (m_spectral[i] != NULL) delete m_spectral[i];

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
 * Updates the precomputation cache for Monte Carlo simulations. In case that
 * the energy boundaries have changed or at least one of the model parameters
 * has changed the method computes a spectral node function which has 100
 * nodes per decade containing the multiplicative spectral model values and
 * stores that into a Monte Carlo cache. This cache is then used by the mc()
 * method for simulations.
 ***************************************************************************/
void GModelSpectralMultiplicative::update_mc_cache(const GEnergy& emin,
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
        GEnergies energies = GEnergies(nodes, m_mc_emin, m_mc_emax, "LOG");

        // Append nodes to spectral function
        for (int i = 0; i < energies.size(); ++i) {
            m_mc_spectrum.append(energies[i], eval(energies[i]));
        }

    } // endif: emin and emax have changed

    // Return
    return;
}
