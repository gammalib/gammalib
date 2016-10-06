/***************************************************************************
 *     GModelSpatialComposite.cpp - Spatial point source model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Domenico Tiziani                                 *
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
 * @file GModelSpatialComposite.cpp
 * @brief Spatial composite model class implementation
 * @author Domenico Tiziani
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GModelSpatialComposite.hpp"
#include "GModelSpatialRegistry.hpp"

/* __ Constants __________________________________________________________ */


/* __ Globals ____________________________________________________________ */
const GModelSpatialComposite g_spatial_comp_seed;
const GModelSpatialRegistry    g_spatial_comp_registry(&g_spatial_comp_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ                 "GModelSpatialComposite::read(GXmlElement&)"
#define G_WRITE                "GModelSpatialComposite::write(GXmlElement&)"
#define G_COMPONENT_INDEX      "GModelSpatialComposite::component(const int&)"
#define G_COMPONENT_NAME       "GModelSpatialComposite::component(const std::string&)"
#define G_APPEND               "GModelSpatialComposite::append(const GModelSpatial&, const std::string&)"

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
 *
 * Constructs empty point source model.
 ***************************************************************************/
GModelSpatialComposite::GModelSpatialComposite(void) : GModelSpatial()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Model type constructor
 *
 * @param[in] dummy Dummy flag.
 * @param[in] type Model type.
 *
 * Constructs empty point source model by specifying a model @p type.
 ***************************************************************************/
GModelSpatialComposite::GModelSpatialComposite(const bool&        dummy,
                                                   const std::string& type) :
                          GModelSpatial()
{
    // Initialise members
    init_members();

    // Set model type
    m_type = type;

    // Return
    return;
}

/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Construct a point source spatial model by extracting information from an
 * XML element. See the read() method for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelSpatialComposite::GModelSpatialComposite(const GXmlElement& xml) :
                          GModelSpatial()
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
 * @param[in] model Spatial composite model.
 ***************************************************************************/
GModelSpatialComposite::GModelSpatialComposite(const GModelSpatialComposite& model) :
    GModelSpatial(model)
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
GModelSpatialComposite::~GModelSpatialComposite(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] model Point source spatial model.
 * @return Point source spatial model.
 ***************************************************************************/
GModelSpatialComposite& GModelSpatialComposite::operator=(const GModelSpatialComposite& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelSpatial::operator=(model);

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
 =                            Public methods                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear point source model
 ***************************************************************************/
void GModelSpatialComposite::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelSpatial::free_members();

    // Initialise members
    this->GModelSpatial::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone point source model
 *
 * @return Pointer to deep copy of point source model.
 ***************************************************************************/
GModelSpatialComposite* GModelSpatialComposite::clone(void) const
{
    // Clone point source model
    return new GModelSpatialComposite(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] photon Incident photon.
 * @param[in] gradients Compute gradients?
 * @return Model value.
 *
 * Evaluates the spatial part for a point source model. It implements a delta
 * function with respect to the coordinates of the source. For numerical
 * reasons, a certain tolerance is accepted (typically 0.1 arcsec, i.e. well
 * below the angular resolution of gamma-ray telescopes).
 *
 * The method will not compute analytical parameter gradients, even if the
 * @p gradients argument is set to true. Point source parameter gradients
 * need to be computed numerically.
 ***************************************************************************/
double GModelSpatialComposite::eval(const GPhoton& photon,
                                      const bool&    gradients) const
{
    // Initalise value
    double value = 0.;
    
    // Sum over all components
    for (unsigned int i = 0; i < components(); i++) {
        value += m_components[i]->eval(photon, gradients);
    }
    
    // Normalise
    value /= m_components.size();
    

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Returns MC sky direction
 *
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @param[in,out] ran Random number generator.
 * @return Sky direction.
 *
 * Draws an arbitrary sky direction for the point source model. As the point
 * source is a point in the sky, the methods always returns the directon of
 * the point source.
 ***************************************************************************/
GSkyDir GModelSpatialComposite::mc(const GEnergy& energy,
                                     const GTime&   time,
                                     GRan&          ran) const
{
    // Choose one component
    
    int idx = 1;
    double u = ran.uniform() * m_components.size();
    
    while (idx < m_components.size()) {
            if (u <= idx) {
                break;
            }
            idx++;
        }
        GModelSpatial *component = m_components[idx-1];

        // Simulate component
        GSkyDir sky_dir = component->mc(energy, time, ran);
        std::cout << sky_dir << std::endl;
    
    // Return sky direction
    return (sky_dir);
}


/***********************************************************************//**
 * @brief Checks where model contains specified sky direction
 *
 * @param[in] dir Sky direction.
 * @param[in] margin Margin to be added to sky direction (degrees)
 * @return True if the model contains the sky direction.
 *
 * Signals whether a sky direction is contained in the point source
 * model.
 ***************************************************************************/
bool GModelSpatialComposite::contains(const GSkyDir& dir,
                                        const double&  margin) const
{
    // Loop over all components
    for (unsigned int idx_comp = 0; idx_comp < components(); idx_comp++) {
        if (m_components[idx_comp]->contains(dir, margin))
            return true;
    }
    return false;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * Reads the spatial information from an XML element.
 ***************************************************************************/
void GModelSpatialComposite::read(const GXmlElement& xml)
{   
    // Get number of spatial components
    int n_comp = xml.elements("spatialModel");
    
    // Loop over spatial elements
    for(int i = 0; i < n_comp; ++i) {
        
        // Get spatial XML element
        const GXmlElement* spec = xml.element("spatialModel", i);
        
        // Initialise a spatial registry object
        GModelSpatialRegistry registry;
        
        // Read spatial model
        GModelSpatial* ptr = registry.alloc(*spec);
        
        // Get component attribute from XML file
        std::string component_name = spec->attribute("component");
        
        // Append spatial component to container
        append(*ptr, component_name);
        
    } // endfor: loop over components
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * @exception GException::model_invalid_spatial
 *            Existing XML element is not of type 'SkyDirFunction'
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Write the point source information into an XML element with the following
 * format
 *
 *     <spatialModel type="PointSource">
 *       <parameter free="0" max="360" min="-360" name="RA" scale="1" value="83.6331" />
 *       <parameter free="0" max="90" min="-90" name="DEC" scale="1" value="22.0145" />
 *     </spatialModel>
 *
 * @todo The case that an existing spatial XML element with "GLON" and "GLAT"
 *       as coordinates is not supported.
 ***************************************************************************/
void GModelSpatialComposite::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "") {
        xml.attribute("type", type());
    }

    // Verify model type
    if (xml.attribute("type") != type()) {
        throw GException::model_invalid_spatial(G_WRITE, xml.attribute("type"),
              "Spatial model is not of type \""+type()+"\".");
    }

    // Verify that XML element is composed from spatial models 
    if (xml.elements("spatialModel") != xml.elements()) {
        throw GException::model_invalid_spatial(G_WRITE, type(),
        "Composite spatial model must be composed solely from spatial models.");
    }

    // Write all model components
    for (unsigned int idx_comp = 0; idx_comp < m_components.size(); idx_comp++) {
        GModelSpatial *component = m_components[idx_comp];
        
        // Find xml element with matching name
        GXmlElement *matching_model = 0;
        for (unsigned int idx_elt = 0; idx_elt < xml.elements("spatialModel"); idx_elt++) {
            if (xml.element("spatialModel", idx_elt)->attribute("name") == 
                m_component_names[idx_comp]) {
                matching_model = xml.element("spatialModel", idx_elt);
                break;
            }
        } // endfor: loop over all xml elements

        // Strip prefix from parameter names

        // Loop over all parameters of component
        for (unsigned int idx_par = 0; idx_par < component->size(); idx_par++) {
            
            // Get parameter
            GModelPar& par = (*component)[idx_par];
            
            // Get name components
            std::vector<std::string> name_components = gammalib::split(par.name(), ":");

            // Remove first component
            name_components.erase(name_components.begin());

            // Join new name
            std::string new_name = name_components[0];
            for (unsigned int i = 1; i < name_components.size(); i++) {
                new_name.append(":");
                new_name.append(name_components[i]);
            }
            
            // Set new name
            par.name(new_name);

        } // endfor looped over all parameters

        // Write component
        if (matching_model) {
            // Use existing xml element
            component->write(*matching_model);
        } else {
            // Create new xml element
            GXmlElement new_element("spatialModel");
            component->write(new_element);
            xml.append(new_element);
        }
        
        // Re-add prefix to parameter names
        
        std::string prefix = m_component_names[idx_comp];

        // Loop over all parameters of component
        for (unsigned int idx_par = 0; idx_par < component->size(); idx_par++) {
            
            // Get parameter
            GModelPar& par = (*component)[idx_par];
            
            // Prepend prefix
            std::string new_name = prefix;
            new_name.append(":");
            new_name.append(par.name());
            
            // Set new name
            par.name(new_name);

        } // endfor looped over all parameters

        
        
    } // endfor: loop over all components

    // TODO: Remove unused xml elements

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return number of model components
 *
 * @return Number of model components
 ***************************************************************************/
int GModelSpatialComposite::components(void) const {
    return m_components.size();
}


/***********************************************************************//**
 * @brief Append spatial component
 *
 * @param[in] model Spatial model component to append
 * @param[in] name Name of spatial model (can be empty)
 * 
 * Appends a spatial component to the composite model
 ***************************************************************************/
void GModelSpatialComposite::append(const GModelSpatial& component,
                                    const std::string& name) {

// Append model container
    m_components.push_back(component.clone());

    // Get index of latest model
    int index = m_components.size()-1;
    
    // Use model index if component name is empty
    std::string component_name = !name.empty() ? name : gammalib::str(m_components.size());
    
    // Check if component name is unique, throw exception if not
    if (gammalib::contains(m_component_names, component_name)) {
    	std::string msg =
            "Attempt to append component with name \""+component_name+"\" to composite"
            " spectral model container, but a component with the same name exists already.\n"
            "Every component in the container needs a unique name. On default the system will"
            "increment an integer if no component name is provided.";
        throw GException::invalid_value(G_APPEND, msg);
    }

    // Add component name (for now simple number)
    m_component_names.push_back(component_name);
    
    // Get number of spectral parameters from model
    int npars = m_components[index]->size();
    
    // Loop over model parameters
    for (int ipar = 0; ipar < npars; ++ipar) {
        
        // Get model parameter
        GModelPar* par = &(m_components[index]->operator[](ipar));
        
        // Modify parameter name
        par->name(component_name+":"+par->name());
        
        // Append model parameter with new name to internal container
        m_pars.push_back(par);
        
    } // endfor: loop over model parameters
}

std::string GModelSpatialComposite::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

    	// Append header
        result.append("=== GModelSpatialComposite ===");

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



/***********************************************************************//**
 * @brief Returns spatial component element
 *
 * @param[in] index Index of spatial component.
 * @return Spatial model.
 *
 * Returns a spatial component to the composite model
 ***************************************************************************/
const GModelSpatial* GModelSpatialComposite::component(const int& index) const
{
	// Check if index is in validity range
	if (index >= m_components.size() || index < 0) {
		throw GException::out_of_range(G_COMPONENT_INDEX, "Component Index", index, m_components.size(),"");
	}

	// Return spectral component
	return m_components[index];

}


/***********************************************************************//**
 * @brief Returns pointer to specific spatial component
 *
 * @param[in] name Name of spatial component.
 * @return Spatial model.
 *
 * Returns a spatial component of the composite model
 ***************************************************************************/
const GModelSpatial* GModelSpatialComposite::component(const std::string& name) const
{
	// Check if model name is found
	int index = -1;
	for(int i = 0; i < m_component_names.size(); ++i) {
            if (m_component_names[i] == name) {
			index = i;
			break;
		}
	}

	// Check if component name was found
	if (index == -1) {
		throw GException::model_not_found(G_COMPONENT_NAME, name,"");
	}

	// Return spectral component
	return m_components[index];

}

/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModelSpatialComposite::init_members(void)
{
    // Initialise model type
    m_type = "SpatialComposite";

    // Initialise models vector
    m_components.clear();
    m_component_names.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Spatial composite model.
 ***************************************************************************/
void GModelSpatialComposite::copy_members(const GModelSpatialComposite& model)
{
    // Copy members
    m_type = model.m_type;
    m_component_names = model.m_component_names;
    m_components.clear();
    m_pars.clear();
    for (unsigned int idx_comp = 0; idx_comp < model.m_components.size(); idx_comp ++) {
        m_components.push_back(model.m_components[idx_comp]->clone());

        // Get number of parameters to append
        int npars = m_components[idx_comp]->size();
    
        // Append parameters
        for (unsigned int idx_par = 0; idx_par < npars; idx_par++) {
            // Get model parameter
            GModelPar& par = (*m_components[idx_comp])[idx_par];
        
            // Append to paramaters
            m_pars.push_back(&par);
        }
    }
    
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatialComposite::free_members(void)
{
    // Free model components
    for (unsigned int i = 0; i < m_components.size(); i ++) {
        
        // Delete component i
        delete m_components[i];

        // Signal free pointer
        m_components[i] = NULL;
    }
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns an unused name for a model component
 * 
 * @return Model name
 ***************************************************************************/
std::string GModelSpatialComposite::unique_component_name(void) {
    unsigned int index = components();
    while (1) {
        // Generate name from name index
        std::string s = gammalib::str(index++);

        // Check if name is unique
        if (component_name_is_unique(s)) { 
            // Return
            return s;
        }
    }
}

/***********************************************************************//**
 * @brief Checks if @p name is not used in component names
 * 
 * @param[in] name Name to check
 * 
 * @return name unique?
 ***************************************************************************/
bool GModelSpatialComposite::component_name_is_unique(std::string name) const {

    // Check if name is unique
    bool is_unique = true;
    for (unsigned int i = 0; i < m_component_names.size(); i++) {
        if (name == m_component_names[i]) {
            is_unique = false;
            break;
        }
    }
    return is_unique;
}
