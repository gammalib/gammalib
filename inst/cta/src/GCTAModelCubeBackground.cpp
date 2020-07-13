/***************************************************************************
 *       GCTAModelCubeBackground.cpp - CTA cube background model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2020 by Michael Mayer                               *
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
 * @file GCTAModelCubeBackground.cpp
 * @brief CTA cube background model class implementation
 * @author Michael Mayer
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
#include "GModelSpectralNodes.hpp"
#include "GCTAModelCubeBackground.hpp"
#include "GCTAObservation.hpp"
#include "GCTAResponseCube.hpp"
#include "GCTACubeBackground.hpp"
#include "GCTASupport.hpp"

/* __ Globals ____________________________________________________________ */
const GCTAModelCubeBackground g_cta_inst_background_seed;
const GModelRegistry          g_cta_inst_background_registry(&g_cta_inst_background_seed);

/* __ Method name definitions ____________________________________________ */
#define G_EVAL "GCTAModelCubeBackground::eval(GEvent&, GObservation&, bool&)"
#define G_NPRED            "GCTAModelCubeBackground::npred(GEnergy&, GTime&,"\
                                                            " GObservation&)"
#define G_MC              "GCTAModelCubeBackground::mc(GObservation&, GRan&)"
#define G_XML_SPECTRAL  "GCTAModelCubeBackground::xml_spectral(GXmlElement&)"
#define G_XML_TEMPORAL  "GCTAModelCubeBackground::xml_temporal(GXmlElement&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_USE_NPRED_CACHE

/* __ Debug definitions __________________________________________________ */
//#define G_DUMP_MC
//#define G_DEBUG_NPRED

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTAModelCubeBackground::GCTAModelCubeBackground(void) : GModelData()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs a CTA cube background model from the information provided by
 * an XML elements (see GCTAModelCubeBackground::read method).
 ***************************************************************************/
GCTAModelCubeBackground::GCTAModelCubeBackground(const GXmlElement& xml) :
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
 * Constructs a CTA cube background model from a @p spectral model component.
 * The temporal component is assumed to be constant.
 ***************************************************************************/
GCTAModelCubeBackground::GCTAModelCubeBackground(const GModelSpectral& spectral) :
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
 * @brief Construct from model components
 *
 * @param[in] spectral Spectral model component.
 * @param[in] temporal Temporal model component.
 *
 * Constructs a CTA cube background model from a @p spectral and a
 * @p temporal component.
 ***************************************************************************/
GCTAModelCubeBackground::GCTAModelCubeBackground(const GModelSpectral& spectral,
                                                 const GModelTemporal& temporal) :
                         GModelData()
{
    // Initialise members
    init_members();

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
 * @param[in] bgd CTA cube background model.
 ***************************************************************************/
GCTAModelCubeBackground::GCTAModelCubeBackground(const GCTAModelCubeBackground& bgd) :
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
GCTAModelCubeBackground::~GCTAModelCubeBackground(void)
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
 * @param[in] bgd CTA cube background model.
 * @return CTA cube background model.
 ***************************************************************************/
GCTAModelCubeBackground& GCTAModelCubeBackground::operator=(const GCTAModelCubeBackground& bgd)
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
 * @brief Clear CTA cube background model
 *
 * This method properly resets the CTA cube background model to an
 * initial state.
 ***************************************************************************/
void GCTAModelCubeBackground::clear(void)
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
 * @brief Clone CTA cube background model
 *
 * @return Pointer to deep copy of CTA cube background model.
 ***************************************************************************/
GCTAModelCubeBackground* GCTAModelCubeBackground::clone(void) const
{
    return new GCTAModelCubeBackground(*this);
}


/***********************************************************************//**
 * @brief Set spectral model component
 *
 * @param[in] spectral Pointer to spectral model component.
 *
 * Sets the spectral model component of the model.
 ***************************************************************************/
void GCTAModelCubeBackground::spectral(const GModelSpectral* spectral)
{
    // Free spectral model component only if it differs from current
    // component. This prevents unintential deallocation of the argument
    if ((m_spectral != NULL) && (m_spectral != spectral)) {
        delete m_spectral;
    }

    // Clone spectral model component if it exists, otherwise set pointer
    // to NULL
    m_spectral = (spectral != NULL) ? spectral->clone() : NULL;

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set temporal model component
 *
 * @param[in] temporal Pointer to temporal model component.
 *
 * Sets the temporal model component of the model.
 ***************************************************************************/
void GCTAModelCubeBackground::temporal(const GModelTemporal* temporal)
{
    // Free temporal model component only if it differs from current
    // component. This prevents unintential deallocation of the argument
    if ((m_temporal != NULL) && (m_temporal != temporal)) {
        delete m_temporal;
    }

    // Clone temporal model component if it exists, otherwise set pointer
    // to NULL
    m_temporal = (temporal != NULL) ? temporal->clone() : NULL;

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return background rate in units of events MeV\f$^{-1}\f$
 *        s\f$^{-1}\f$ sr\f$^{-1}\f$
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 * @param[in] gradients Compute gradients?
 * @return Background rate (events MeV\f$^{-1}\f$ s\f$^{-1}\f$ sr\f$^{-1}\f$).
 *
 * Evaluates the background model. The method returns a real rate, defined
 * as the number of counts per MeV, steradian and ontime.
 *
 * If the @p gradients flag is true the method will also set the parameter
 * gradients of the model parameters.
 ***************************************************************************/
double GCTAModelCubeBackground::eval(const GEvent&       event,
                                     const GObservation& obs,
                                     const bool&         gradients) const
{
    // Get reference on CTA response cube from observation and reference on
    // CTA instrument direction from event
    const GCTAResponseCube& rsp = gammalib::cta_rsp_cube(G_EVAL, obs);
    const GCTAInstDir&      dir = gammalib::cta_dir(G_EVAL, event);

    // Retrieve reference to CTA cube background
    const GCTACubeBackground& bgd = rsp.background();

    // Evaluate function
    double spat = bgd(dir, event.energy());
    double spec = (spectral() != NULL)
                  ? spectral()->eval(event.energy(), event.time(), gradients)
                  : 1.0;
    double temp = (temporal() != NULL)
                  ? temporal()->eval(event.time(), gradients) : 1.0;

    // Compute value. Note that background rates are already per
    // livetime, hence no deadtime correction is needed here.
    double value = spat * spec * temp;

    // If gradients were requested then multiply factors to spectral and
    // temporal gradients
    if (gradients) {

        // Multiply factors to spectral gradients
        if (spectral() != NULL) {
            double fact = spat * temp;
            if (fact != 1.0) {
                for (int i = 0; i < spectral()->size(); ++i)
                    (*spectral())[i].factor_gradient((*spectral())[i].factor_gradient() * fact);
            }
        }

        // Multiply factors to temporal gradients
        if (temporal() != NULL) {
            double fact = spat * spec;
            if (fact != 1.0) {
                for (int i = 0; i < temporal()->size(); ++i)
                    (*temporal())[i].factor_gradient((*temporal())[i].factor_gradient() * fact);
            }
        }

    } // endif: gradients were requested

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Return spatially integrated background rate in units of
 *        events MeV\f$^{-1}\f$ s\f$^{-1}\f$
 *
 * @param[in] obsEng Measured event energy.
 * @param[in] obsTime Measured event time.
 * @param[in] obs Observation.
 * @return Spatially integrated background rate
 *         (events MeV\f$^{-1}\f$ s\f$^{-1}\f$)
 *
 * Spatially integrates the cube background model for a given measured event
 * energy and event time. The method returns a real rate,
 * defined as the number of counts per MeV and ontime.
 ***************************************************************************/
double GCTAModelCubeBackground::npred(const GEnergy&      obsEng,
                                      const GTime&        obsTime,
                                      const GObservation& obs) const
{
    // Initialise result
    double npred     = 0.0;
    bool   has_npred = false;

    // Build unique identifier
    std::string id = obs.instrument() + "::" + obs.id();

    // Check if Npred value is already in cache
    #if defined(G_USE_NPRED_CACHE)
    if (!m_npred_names.empty()) {

        // Search for unique identifier, and if found, recover Npred value
        // and break
        for (int i = 0; i < m_npred_names.size(); ++i) {
            if (m_npred_names[i] == id && m_npred_energies[i] == obsEng) {
                npred     = m_npred_values[i];
                has_npred = true;
                #if defined(G_DEBUG_NPRED)
                std::cout << "GCTAModelCubeBackground::npred:";
                std::cout << " cache=" << i;
                std::cout << " npred=" << npred << std::endl;
                #endif
                break;
            }
        }

    } // endif: there were values in the Npred cache
    #endif

    // Continue only if no Npred cache value has been found
    if (!has_npred) {

        // Evaluate only if model is valid
        if (valid_model()) {

            // Get reference on CTA response cube from observation
            const GCTAResponseCube& rsp = gammalib::cta_rsp_cube(G_NPRED, obs);

            // Retrieve CTA background
            const GCTACubeBackground& bgd = rsp.background();

            // Get log10 of energy in TeV
            double logE = obsEng.log10TeV();

            // Integrate the background map at a certain energy
            npred = bgd.integral(logE);

            // Store result in Npred cache
            #if defined(G_USE_NPRED_CACHE)
            m_npred_names.push_back(id);
            m_npred_energies.push_back(obsEng);
            m_npred_times.push_back(obsTime);
            m_npred_values.push_back(npred);
            #endif

            // Debug: Check for NaN
            #if defined(G_NAN_CHECK)
            if (gammalib::is_notanumber(npred) || gammalib::is_infinite(npred)) {
                std::string origin  = "GCTAModelCubeBackground::npred";
                std::string message = " NaN/Inf encountered (npred=" +
                                      gammalib::str(npred) + ")";
                gammalib::warning(origin, message);
            }
            #endif

        } // endif: model was valid

    } // endif: Npred computation required

    // Multiply in spectral and temporal components
    npred *= spectral()->eval(obsEng, obsTime);
    npred *= temporal()->eval(obsTime);

    // Return Npred
    return npred;
}


/***********************************************************************//**
 * @brief Return simulated list of events
 *
 * @param[in] obs Observation.
 * @param[in] ran Random number generator.
 * @return Pointer to list of simulated events (needs to be de-allocated by
 *         client)
 *
 * @exception GException::feature_not_implemented
 *            Specified observation is not a CTA observation.
 *
 * The simulation of an event list from a cube background model is not
 * implemented, hence the method will always throw an exception.
 ***************************************************************************/
GCTAEventList* GCTAModelCubeBackground::mc(const GObservation& obs, GRan& ran) const
{
    // Feature not yet implemented
    throw GException::feature_not_implemented(G_MC,
          "MC computation not implemented for binned analysis.");

    // Return NULL pointer
    return NULL;

}


/***********************************************************************//**
 * @brief Read CTA cube background model from XML element
 *
 * @param[in] xml XML element.
 *
 * Set up CTA cube background model from the information provided by
 * an XML element. The XML element is expected to have the following
 * structure
 *
 *     <source name="..." type="..." instrument="...">
 *       <spectrum type="...">
 *         ...
 *       </spectrum>
 *     </source>
 *
 * Optionally, a temporal model may be provided using the following
 * syntax
 *
 *     <source name="..." type="..." instrument="...">
 *       <spectrum type="...">
 *         ...
 *       </spectrum>
 *       <temporalModel type="...">
 *         ...
 *       </temporalModel>
 *     </source>
 *
 * If no temporal component is found a constant model is assumed.
 ***************************************************************************/
void GCTAModelCubeBackground::read(const GXmlElement& xml)
{
    // Clear model
    clear();

    // Initialise XML elements
    const GXmlElement* spectral = NULL;

    // Get pointer on spectrum
    spectral = xml.element("spectrum", 0);

    // Extract spectral model
    m_spectral = xml_spectral(*spectral);

    // Optionally get temporal model
    if (xml.elements("temporalModel")) {
        const GXmlElement* temporal = xml.element("temporalModel", 0);
        m_temporal = xml_temporal(*temporal);
    }
    else {
        GModelTemporalConst constant;
        m_temporal = constant.clone();
    }

    // Read model attributes
    read_attributes(xml);

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write CTA cube background model into XML element
 *
 * @param[in] xml XML element.
 *
 * Write CTA cube background model information into an XML element.
 * The XML element will have the following structure
 *
 *     <source name="..." type="..." instrument="...">
 *       <spectrum type="...">
 *         ...
 *       </spectrum>
 *     </source>
 *
 * If the model contains a non-constant temporal model, the temporal
 * component will also be written following the syntax
 *
 *     <source name="..." type="..." instrument="...">
 *       <spectrum type="...">
 *         ...
 *       </spectrum>
 *       <temporalModel type="...">
 *         ...
 *       </temporalModel>
 *     </source>
 *
 * If no temporal component is found a constant model is assumed.
 ***************************************************************************/
void GCTAModelCubeBackground::write(GXmlElement& xml) const
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

    // If we have a temporal model that is either not a constant, or a
    // constant with a normalization value that differs from 1.0 then
    // write the temporal component into the XML element. This logic
    // assures compatibility with the Fermi/LAT format as this format
    // does not handle temporal components.
    bool write_temporal = ((m_temporal != NULL) &&
                           (m_temporal->type() != "Constant" ||
                            (*m_temporal)[0].value() != 1.0));

    // If no source with corresponding name was found then append one
    if (src == NULL) {
        src = xml.append("source");
        if (spectral() != NULL) src->append(GXmlElement("spectrum"));
        if (write_temporal)     src->append(GXmlElement("temporalModel"));
    }

    // Write spectral model
    if (spectral() != NULL) {
        GXmlElement* spec = src->element("spectrum", 0);
        spectral()->write(*spec);
    }

    // Optionally write temporal model
    if (write_temporal) {
        if (dynamic_cast<GModelTemporalConst*>(temporal()) == NULL) {
            GXmlElement* temp = src->element("temporalModel", 0);
            temporal()->write(*temp);
        }
    }

    // Write model attributes
    write_attributes(*src);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print CTA cube background model information
 *
 * @param[in] chatter Chattiness.
 * @return String containing CTA cube background model information.
 ***************************************************************************/
std::string GCTAModelCubeBackground::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAModelCubeBackground ===");

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
void GCTAModelCubeBackground::init_members(void)
{
    // Initialise members
    m_spectral = NULL;
    m_temporal = NULL;

    // Initialise Npred cache
    m_npred_names.clear();
    m_npred_energies.clear();
    m_npred_times.clear();
    m_npred_values.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] bgd CTA background model.
 ***************************************************************************/
void GCTAModelCubeBackground::copy_members(const GCTAModelCubeBackground& bgd)
{
    // Copy cache
    m_npred_names    = bgd.m_npred_names;
    m_npred_energies = bgd.m_npred_energies;
    m_npred_times    = bgd.m_npred_times;
    m_npred_values   = bgd.m_npred_values;

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
void GCTAModelCubeBackground::free_members(void)
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
void GCTAModelCubeBackground::set_pointers(void)
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
bool GCTAModelCubeBackground::valid_model(void) const
{
    // Set result
    bool result = ((spectral() != NULL) && (temporal() != NULL));

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Return pointer to spectral model from XML element
 *
 * @param[in] spectral XML element.
 * @return Pointer to spectral model.
 *
 * Returns pointer to spectral model that is defined in an XML element.
 ***************************************************************************/
GModelSpectral* GCTAModelCubeBackground::xml_spectral(const GXmlElement& spectral) const
{
    // Get spectral model
    GModelSpectralRegistry registry;
    GModelSpectral*        ptr = registry.alloc(spectral);

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Return pointer to temporal model from XML element
 *
 * @param[in] temporal XML element.
 * @return Pointer to temporal model.
 *
 * Returns pointer to temporal model that is defined in an XML element.
 ***************************************************************************/
GModelTemporal* GCTAModelCubeBackground::xml_temporal(const GXmlElement& temporal) const
{
    // Get temporal model
    GModelTemporalRegistry registry;
    GModelTemporal*        ptr = registry.alloc(temporal);

    // Return pointer
    return ptr;
}

