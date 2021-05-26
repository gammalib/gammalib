/***************************************************************************
 * GCTAModelSpatialMultiplicative.cpp - Multiplicative spatial model class *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2018 by Jurgen Knodlseder                                *
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
 * @file GCTAModelSpatialMultiplicative.cpp
 * @brief Multiplicative spatial model class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GRan.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GXmlElement.hpp"
#include "GCTAObservation.hpp"
#include "GCTAInstDir.hpp"
#include "GCTARoi.hpp"
#include "GCTAModelSpatialMultiplicative.hpp"
#include "GCTAModelSpatialRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GCTAModelSpatialMultiplicative g_cta_spatial_multi_seed;
const GCTAModelSpatialRegistry       g_cta_spatial_multi_registry(&g_cta_spatial_multi_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ           "GCTAModelSpatialMultiplicative::read(GXmlElement&)"
#define G_WRITE         "GCTAModelSpatialMultiplicative::write(GXmlElement&)"
#define G_COMPONENT_INDEX   "GCTAModelSpatialMultiplicative::component(int&)"
#define G_COMPONENT_NAME          "GCTAModelSpatialMultiplicative::component"\
                                                             "(std::string&)"
#define G_APPEND "GCTAModelSpatialMultiplicative::append(GCTAModelSpatial&, "\
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
GCTAModelSpatialMultiplicative::GCTAModelSpatialMultiplicative(void) :
                                GCTAModelSpatial()
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
 * Creates instance of a multiplicative spatial model by extracting
 * information from an XML element.
 * See GCTAModelSpatialMultiplicative::read() for more information about the
 * expected structure of the XML element.
 ***************************************************************************/
GCTAModelSpatialMultiplicative::GCTAModelSpatialMultiplicative(const GXmlElement& xml) :
                                GCTAModelSpatial()
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
 * @param[in] model Spatial model.
 ***************************************************************************/
GCTAModelSpatialMultiplicative::GCTAModelSpatialMultiplicative(const GCTAModelSpatialMultiplicative& model) :
                                GCTAModelSpatial()
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
GCTAModelSpatialMultiplicative::~GCTAModelSpatialMultiplicative(void)
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
 * @param[in] model Spatial model.
 ***************************************************************************/
GCTAModelSpatialMultiplicative& GCTAModelSpatialMultiplicative::operator=(const GCTAModelSpatialMultiplicative& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GCTAModelSpatial::operator=(model);

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
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 ***************************************************************************/
void GCTAModelSpatialMultiplicative::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GCTAModelSpatial::free_members();

    // Initialise members
    this->GCTAModelSpatial::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 ***************************************************************************/
GCTAModelSpatialMultiplicative* GCTAModelSpatialMultiplicative::clone(void) const
{
    return new GCTAModelSpatialMultiplicative(*this);
}


/***********************************************************************//**
 * @brief Evaluate function (in units of sr^-1)
 *
 * @param[in] dir Event direction.
 * @param[in] energy Event energy.
 * @param[in] time Event time.
 * @param[in] gradients Compute gradients?
 * @return Model value
 *
 * Evaluate spatial model for a given event direction.
 ***************************************************************************/
double GCTAModelSpatialMultiplicative::eval(const GCTAInstDir& dir,
                                            const GEnergy&     energy,
                                            const GTime&       time,
                                            const bool&        gradients) const
{
    // Initialise result
    double value = 0.0;

    // Check for available spatial components
    if (m_spatial.size() > 0) {

        // Initialise vector of values
        std::vector<double> values(m_spatial.size(), 0.0);

        // Compute values of components and function value
        for (int i = 0; i < m_spatial.size(); ++i) {
            values[i] = m_spatial[i]->eval(dir, energy, time, gradients);
            if (i == 0) {
                value = values[0];
            }
            else {
                value *= values[i];
            }
        }

        // Modify gradients if requested
        if (gradients) {

            // Loop over model components
            for (int i = 0; i < m_spatial.size(); ++i) {

                // Initialise scaling factor
                double factor = 1.0;

                // Loop over other model components and compute factor
                for (int j = 0; j < m_spatial.size(); ++j) {
                    if (i != j) {
                        factor *= values[j];
                    }
                }

                // Loop over model parameters
                for (int ipar = 0; ipar < m_spatial[i]->size(); ++ipar) {

                    // Get reference to model parameter
                    GModelPar& par = m_spatial[i]->operator[](ipar);

                    // Scale parameter gradient
                    par.gradient(par.gradient()*factor);

                } // endfor: loop over model parameters

            } // endfor: loop over models

        } //endif: compute gradients

    } // endif: there were spatial components

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GCTAModelSpatialMultiplicative::eval";
        std::cout << "(dir=" << dir.print();
        std::cout << ", energy=" << energy;
        std::cout << ", time=" << time << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Return maximum function value for Monte Carlo simulations
 *
 * @param[in] obs CTA Observation.
 * @return Maximum function value for Monte Carlo simulations.
 ***************************************************************************/
double GCTAModelSpatialMultiplicative::mc_max_value(const GCTAObservation& obs) const
{
    // Initialise maximum value
    double value = 1.0;

    // Compute maximum value from components
    for (int i = 0; i < m_spatial.size(); ++i) {
        value *= m_spatial[i]->mc_max_value(obs);
    }

    // Return maximum value
    return value;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * Reads the spatial information from an XML element.
 ***************************************************************************/
void GCTAModelSpatialMultiplicative::read(const GXmlElement& xml)
{
    // Get number of spatial components
    int n_spatials = xml.elements("spatialModel");

    // Loop over spatial elements
    for (int i = 0; i < n_spatials; ++i) {

        // Get spatial XML element
        const GXmlElement* spatial = xml.element("spatialModel", i);

        // Allocate a spatial registry object
        GCTAModelSpatialRegistry registry;

        // Read spatial model
        GCTAModelSpatial* ptr = registry.alloc(*spatial);

        // Get component attribute from XML file
        std::string component_name = spatial->attribute("component");

        // Append spatial component to container
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
 * Writes the spatial information into an XML element.
 ***************************************************************************/
void GCTAModelSpatialMultiplicative::write(GXmlElement& xml) const
{
    // Check model type
    gammalib::xml_check_type(G_WRITE, xml, type());

    // Loop over model components
    for (int i = 0; i < m_spatial.size(); i++) {

        // Write spatial model
        if (m_spatial[i] != NULL) {

            // Create new spectrum node
            xml.append(GXmlElement("spatialModel"));

            // Get new spatial node
            GXmlElement* spatial = xml.element("spatialModel",
                                               xml.elements("spatialModel")-1);

            // Create temporary copy of the spatial model. This is a kluge to
            // write out the original parameters.
            GCTAModelSpatial* cpy = m_spatial[i]->clone();

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

            // Write spatial component
            cpy->write(*spatial);

            // Add component name if previously available
            if (m_components[i] != gammalib::str(i+1)) {
                spatial->attribute("component", m_components[i]);
            }

            // Remove temporary copy
            delete cpy;

        } // endif: spatial model was not NULL

    } // endfor: loop over model components

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append spatial component
 *
 * @param[in] spatial Spatial model component.
 * @param[in] name Name of spatial component (can be empty).
 *
 * @exception GException::invalid_value
 *            Invalid component name specified
 *
 * Appends a spatial component to the multiplicative model
 ***************************************************************************/
void GCTAModelSpatialMultiplicative::append(const GCTAModelSpatial& spatial,
                                            const std::string&      name)
{
    // Append model container
    m_spatial.push_back(spatial.clone());

    // Get index of latest model
    int index = m_spatial.size()-1;

    // Use model index as component name if component name is empty
    std::string component_name = !name.empty() ? name
                                               : gammalib::str(m_spatial.size());

    // Check if component name is unique, throw exception if not
    if (gammalib::contains(m_components, component_name)) {
        std::string msg = "Attempt to append component with name \""+
                          component_name+"\" to multiplicative spatial model "
                          "container, but a component with the same name exists "
                          "already. Every component in the container needs a "
                          "unique name. On default the system will increment "
                          "an integer if no component name is provided.";
        throw GException::invalid_value(G_APPEND, msg);
    }

    // Add component name (for now simple number)
    m_components.push_back(component_name);

    // Get number of spatial parameters from model
    int npars = m_spatial[index]->size();

    // Loop over model parameters
    for (int ipar = 0; ipar < npars; ++ipar) {

        // Get model parameter
        GModelPar* par = &(m_spatial[index]->operator[](ipar));

        // Prepend component name to parameter name
        par->name(component_name+":"+par->name());

        // Append model parameter with new name to internal container
        m_pars.push_back(par);

    } // endfor: loop over model parameters

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return spatial model component by index
 *
 * @param[in] index Index of spectral component.
 * @return Pointer to spatial model.
 *
 * @exception GException::out_of_range
 *            Model index out of valid range
 *
 * Returns a component of the multiplicative spatial model by @p index.
 ***************************************************************************/
const GCTAModelSpatial* GCTAModelSpatialMultiplicative::component(const int& index) const
{
    // Check if index is in validity range
    if (index >= m_spatial.size() || index < 0) {
        throw GException::out_of_range(G_COMPONENT_INDEX, "Component Index",
                                       index, m_spatial.size());
    }

    // Return spatial component
    return m_spatial[index];
}


/***********************************************************************//**
 * @brief Return spatial model component by name
 *
 * @param[in] name Name of spatial component.
 * @return Pointer to spatial model.
 *
 * @exception GException::invalid_argument
 *            Spatial component not found
 *
 * Returns a component of the multiplicative spatial model by @p name.
 ***************************************************************************/
const GCTAModelSpatial* GCTAModelSpatialMultiplicative::component(const std::string& name) const
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
        std::string msg = "Model component \""+name+"\" not found in "
                          "multiplicative spatial model.";
        throw GException::invalid_argument(G_COMPONENT_NAME, msg);
    }

    // Return spatial component
    return m_spatial[index];
}


/***********************************************************************//**
 * @brief Print multiplicative spatial model information
 *
 * @param[in] chatter Chattiness.
 * @return String containing model information.
 ***************************************************************************/
std::string GCTAModelSpatialMultiplicative::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAModelSpatialMultiplicative ===");

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
void GCTAModelSpatialMultiplicative::init_members(void)
{
    // Initialise model type
    m_type = "Multiplicative";

    // Clear spatial models
    m_spatial.clear();
    m_components.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Spatial background model component.
 ***************************************************************************/
void GCTAModelSpatialMultiplicative::copy_members(const GCTAModelSpatialMultiplicative& model)
{
    // Copy members
    m_type       = model.m_type;
    m_components = model.m_components;

    // Copy pointer(s) of spatial component
    m_spatial.clear();
    for (int i = 0; i < model.components(); ++i) {
    	m_spatial.push_back(model.m_spatial[i]->clone());
    }

    // Store pointers to spatial parameters
    m_pars.clear();
    for (int i = 0; i < model.components(); ++i) {

        // Retrieve spatial model
        GCTAModelSpatial* spatial = m_spatial[i];

        // Loop over parameters
        for (int ipar = 0; ipar < spatial->size(); ++ipar) {

            // Get model parameter reference
            GModelPar& par = spatial->operator[](ipar);

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
void GCTAModelSpatialMultiplicative::free_members(void)
{
    // Free memory
    for (int i = 0; i < m_spatial.size(); ++i) {

        // Delete component i
        if (m_spatial[i] != NULL) delete m_spatial[i];

        // Signal free pointer
        m_spatial[i] = NULL;
    
    }

    // Return
    return;
}
