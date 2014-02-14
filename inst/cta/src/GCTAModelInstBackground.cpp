/***************************************************************************
 *   GCTAModelInstBackground.cpp - CTA instrument background model class   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Juergen Knoedlseder                              *
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
 * @file GCTAModelInstBackground.cpp
 * @brief CTA instrument background model class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GModelRegistry.hpp"
#include "GModelSpectralRegistry.hpp"
#include "GModelTemporalRegistry.hpp"
#include "GModelTemporalConst.hpp"
#include "GCTAModelInstBackground.hpp"
#include "GCTAObservation.hpp"

/* __ Globals ____________________________________________________________ */
const GCTAModelInstBackground g_cta_inst_background_seed;
const GModelRegistry          g_cta_inst_background_registry(&g_cta_inst_background_seed);

/* __ Method name definitions ____________________________________________ */
#define G_EVAL        "GCTAModelInstBackground::eval(GEvent&, GObservation&)"
#define G_EVAL_GRADIENTS   "GCTAModelInstBackground::eval_gradients(GEvent&,"\
                                                            " GObservation&)"
#define G_XML_SPECTRAL  "GCTAModelInstBackground::xml_spectral(GXmlElement&)"
#define G_XML_TEMPORAL  "GCTAModelInstBackground::xml_temporal(GXmlElement&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTAModelInstBackground::GCTAModelInstBackground(void) : GModelData()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs a CTA instrumental background model from the information that
 * is found in a XML element. Please refer to the method
 * GCTAModelInstBackground::read to learn more about the information that is
 * expected in the XML element.
 ***************************************************************************/
GCTAModelInstBackground::GCTAModelInstBackground(const GXmlElement& xml) :
                         GModelData(xml)
{
    // Initialise members
    init_members();

    // Read XML
    read(xml);

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Construct from spectral component
 *
 * @param[in] spectral Spectral model component.
 *
 * Constructs a CTA instrumental background model from a spectral
 * model component. The temporal component is assumed to be constant.
 * Please refer to the classe GModelSpectral to learn more about the
 * definition of the spectral components.
 ***************************************************************************/
GCTAModelInstBackground::GCTAModelInstBackground(const GModelSpectral& spectral) :
                         GModelData()
{
    // Initialise members
    init_members();

    // Allocate temporal constant model
    GModelTemporalConst temporal;

    // Clone model components
    m_spectral = spectral.clone();
    m_temporal = temporal.clone();

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] bgd CTA instrument background model.
 ***************************************************************************/
GCTAModelInstBackground::GCTAModelInstBackground(const GCTAModelInstBackground& bgd) :
                         GModelData(bgd)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(bgd);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAModelInstBackground::~GCTAModelInstBackground(void)
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
 * @param[in] bgd CTA instrument background model.
 * @return CTA instrument background model.
 ***************************************************************************/
GCTAModelInstBackground& GCTAModelInstBackground::operator=(const GCTAModelInstBackground& bgd)
{
    // Execute only if object is not identical
    if (this != &bgd) {

        // Copy base class members
        this->GModelData::operator=(bgd);

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(bgd);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear CTA instrument background model
 *
 * This method properly resets the CTA instrument background model to an
 * initial state.
 ***************************************************************************/
void GCTAModelInstBackground::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelData::free_members();

    // Initialise members
    this->GModelData::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone CTA instrument background model
 *
 * @return Pointer to deep copy of CTA instrument background model.
 ***************************************************************************/
GCTAModelInstBackground* GCTAModelInstBackground::clone(void) const
{
    return new GCTAModelInstBackground(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 * @return Function value.
 *
 * @exception GException::invalid_argument
 *            Specified observation is not of the expected type.
 *
 * @todo Implement method
 ***************************************************************************/
double GCTAModelInstBackground::eval(const GEvent& event,
                                     const GObservation& obs) const
{
    // Get pointer on CTA observation
    const GCTAObservation* ctaobs = dynamic_cast<const GCTAObservation*>(&obs);
    if (ctaobs == NULL) {
        std::string msg = "Specified observation is not a CTA observation.\n" +
                          obs.print();
        throw GException::invalid_argument(G_EVAL, msg);
    }

    // Return
    return 0.0;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 * @return Function value.
 *
 * @exception GException::invalid_argument
 *            Specified observation is not of the expected type.
 *
 * @todo Implement method
 ***************************************************************************/
double GCTAModelInstBackground::eval_gradients(const GEvent& event,
                                               const GObservation& obs) const
{
    // Get pointer on CTA observation
    const GCTAObservation* ctaobs = dynamic_cast<const GCTAObservation*>(&obs);
    if (ctaobs == NULL) {
        std::string msg = "Specified observation is not a CTA observation.\n" +
                          obs.print();
        throw GException::invalid_argument(G_EVAL_GRADIENTS, msg);
    }

    // Return
    return 0.0;
}



/***********************************************************************//**
 * @brief Return spatially integrated data model
 *
 * @param[in] obsEng Measured event energy.
 * @param[in] obsTime Measured event time.
 * @param[in] obs Observation.
 * @return Spatially integrated model.
 *
 * @todo Implement method
 ***************************************************************************/
double GCTAModelInstBackground::npred(const GEnergy&      obsEng,
                                      const GTime&        obsTime,
                                      const GObservation& obs) const
{
    // Return
    return 0.0;
}


/***********************************************************************//**
 * @brief Return simulated list of events
 *
 * @param[in] obs Observation.
 * @param[in] ran Random number generator.
 * @return Pointer to list of simulated events (needs to be de-allocated by
 *         client)
 *
 * @todo Implement method
 ***************************************************************************/
GCTAEventList* GCTAModelInstBackground::mc(const GObservation& obs, GRan& ran) const
{
    // Return
    return NULL;
}


/***********************************************************************//**
 * @brief Read CTA instrument background model from XML element
 *
 * @param[in] xml XML element.
 *
 * The model is composed of a spectrum component ('spectrum') and optionally
 * of a temporal component ('temporalModel'). If no temporal component is
 * found a constant model is assumed.
 ***************************************************************************/
void GCTAModelInstBackground::read(const GXmlElement& xml)
{
    // Clear model
    clear();

    // Initialise XML elements
    const GXmlElement* spectral = NULL;
    const GXmlElement* temporal = NULL;

    // Get pointer on spectrum
    spectral = xml.element("spectrum", 0);

    // Extract spectral model
    m_spectral = xml_spectral(*spectral);

    // Optionally get temporal model
    if (xml.elements("temporalModel")) {
        temporal   = xml.element("temporalModel", 0);
        m_temporal = xml_temporal(*temporal);
    }
    else {
        GModelTemporalConst constant;
        m_temporal = constant.clone();
    }

    // Set model name
    name(xml.attribute("name"));

    // Set instruments
    instruments(xml.attribute("instrument"));

    // Set observation identifiers
    ids(xml.attribute("id"));

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write CTA instrument background model into XML element
 *
 * @param[in] xml XML element.
 *
 * @todo Document method.
 ***************************************************************************/
void GCTAModelInstBackground::write(GXmlElement& xml) const
{
    // Initialise pointer on source
    GXmlElement* src = NULL;

    // Search corresponding source
    int n = xml.elements("source");
    for (int k = 0; k < n; ++k) {
        GXmlElement* element = xml.element("source", k);
        if (element->attribute("name") == name()) {
            src = element;
            break;
        }
    }

    // If no source with corresponding name was found then append one
    if (src == NULL) {
        src = xml.append("source");
        if (spectral() != NULL) src->append(GXmlElement("spectrum"));
        //if (temporal() != NULL) src->append(GXmlElement("temporalModel"));
    }

    // Set model type, name and optionally instruments
    src->attribute("name", name());
    src->attribute("type", type());
    if (instruments().length() > 0) {
        src->attribute("instrument", instruments());
    }
    std::string identifiers = ids();
    if (identifiers.length() > 0) {
        src->attribute("id", identifiers);
    }

    // Write spectral model
    if (spectral() != NULL) {
        GXmlElement* spec = src->element("spectrum", 0);
        spectral()->write(*spec);
    }

    // Write temporal model
    /*
    if (temporal() != NULL) {
        if (dynamic_cast<GModelTemporalConst*>(temporal()) == NULL) {
            GXmlElement* temp = src->element("temporalModel", 0);
            temporal()->write(*temp);
        }
    }
    */

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print CTA instrument background model information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing CTA instrument background model information.
 ***************************************************************************/
std::string GCTAModelInstBackground::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAModelInstBackground ===");

        // Determine number of parameters per type
        int n_spectral = (spectral() != NULL) ? spectral()->size() : 0;
        int n_temporal = (temporal() != NULL) ? temporal()->size() : 0;

        // Append attributes
        result.append("\n"+print_attributes());

        // Append model type
        result.append("\n"+gammalib::parformat("Model type"));
        if (n_spectral > 0) {
            result.append("\""+spectral()->type()+"\"");
            if (n_temporal > 0) {
                result.append(" * ");
            }
        }
        if (n_temporal > 0) {
            result.append("\""+temporal()->type()+"\"");
        }

        // Append parameters
        result.append("\n"+gammalib::parformat("Number of parameters") +
                      gammalib::str(size()));
        result.append("\n"+gammalib::parformat("Number of spectral par's") +
                      gammalib::str(n_spectral));
        for (int i = 0; i < n_spectral; ++i) {
            result.append("\n"+(*spectral())[i].print());
        }
        result.append("\n"+gammalib::parformat("Number of temporal par's") +
                      gammalib::str(n_temporal));
        for (int i = 0; i < n_temporal; ++i) {
            result.append("\n"+(*temporal())[i].print());
        }

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCTAModelInstBackground::init_members(void)
{
    // Initialise members
    m_spectral = NULL;
    m_temporal = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] bgd CTA background model.
 ***************************************************************************/
void GCTAModelInstBackground::copy_members(const GCTAModelInstBackground& bgd)
{
    // Clone spectral and temporal model components
    m_spectral = (bgd.m_spectral != NULL) ? bgd.m_spectral->clone() : NULL;
    m_temporal = (bgd.m_temporal != NULL) ? bgd.m_temporal->clone() : NULL;

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAModelInstBackground::free_members(void)
{
    // Free memory
    if (m_spectral != NULL) delete m_spectral;
    if (m_temporal != NULL) delete m_temporal;

    // Signal free pointers
    m_spectral = NULL;
    m_temporal = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set pointers
 *
 * Set pointers to all model parameters. The pointers are stored in a vector
 * that is member of the GModelData base class.
 ***************************************************************************/
void GCTAModelInstBackground::set_pointers(void)
{
    // Clear parameters
    m_pars.clear();

    // Determine the number of parameters
    int n_spectral = (spectral() != NULL) ? spectral()->size() : 0;
    int n_temporal = (temporal() != NULL) ? temporal()->size() : 0;
    int n_pars     = n_spectral + n_temporal;

    // Continue only if there are parameters
    if (n_pars > 0) {

        // Gather spectral parameters
        for (int i = 0; i < n_spectral; ++i) {
            m_pars.push_back(&((*spectral())[i]));
        }

        // Gather temporal parameters
        for (int i = 0; i < n_temporal; ++i) {
            m_pars.push_back(&((*temporal())[i]));
        }

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Verifies if model has all components
 *
 * Returns 'true' if models has a spectral and a temporal component.
 * Otherwise returns 'false'.
 ***************************************************************************/
bool GCTAModelInstBackground::valid_model(void) const
{
    // Set result
    bool result = ((spectral() != NULL) && (temporal() != NULL));

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Construct spectral model from XML element
 *
 * @param[in] spectral XML element containing spectral model information.
 *
 * @exception GException::model_invalid_spectral
 *            Invalid spectral model type encountered.
 *
 * Returns pointer to a spectral model that is defined in an XML element.
 ***************************************************************************/
GModelSpectral* GCTAModelInstBackground::xml_spectral(const GXmlElement& spectral) const
{
    // Get spectral model type
    std::string type = spectral.attribute("type");

    // Get spectral model
    GModelSpectralRegistry registry;
    GModelSpectral*        ptr = registry.alloc(type);

    // If model if valid then read model from XML file
    if (ptr != NULL) {
        ptr->read(spectral);
    }

    // ... otherwise throw an exception
    else {
        throw GException::model_invalid_spectral(G_XML_SPECTRAL, type);
    }

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Construct temporal model from XML element
 *
 * @param[in] temporal XML element containing temporal model information.
 *
 * @exception GException::model_invalid_temporal
 *            Invalid temporal model type encountered.
 *
 * Returns pointer to a temporal model that is defined in an XML element.
 ***************************************************************************/
GModelTemporal* GCTAModelInstBackground::xml_temporal(const GXmlElement& temporal) const
{
    // Get temporal model type
    std::string type = temporal.attribute("type");

    // Get temporal model
    GModelTemporalRegistry registry;
    GModelTemporal*        ptr = registry.alloc(type);

    // If model if valid then read model from XML file
    if (ptr != NULL) {
        ptr->read(temporal);
    }

    // ... otherwise throw an exception
    else {
        throw GException::model_invalid_temporal(G_XML_TEMPORAL, type);
    }

    // Return pointer
    return ptr;
}
