/***************************************************************************
 *                GCOMModelDRM.cpp - COMPTEL DRM model class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2021 by Juergen Knoedlseder                              *
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
 * @file GCOMModelDRM.cpp
 * @brief COMPTEL DRM model class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <typeinfo>
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelRegistry.hpp"
#include "GCOMModelDRM.hpp"
#include "GCOMObservation.hpp"
#include "GCOMEventBin.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GCOMModelDRM   g_com_drm_seed;
const GModelRegistry g_com_drm_registry(&g_com_drm_seed);

/* __ Method name definitions ____________________________________________ */
#define G_EVAL            "GCOMModelDRM::eval(GEvent&, GObservation&, bool&)"
#define G_NRED         "GCOMModelDRM::npred(GEnergy&, GTime&, GObservation&)"
#define G_MC                         "GCOMModelDRM::mc(GObservation&, GRan&)"
#define G_READ                             "GCOMModelDRM::read(GXmlElement&)"
#define G_WRITE                           "GCOMModelDRM::write(GXmlElement&)"

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
 * Constructs an empty COMPTEL DRM model.
 ***************************************************************************/
GCOMModelDRM::GCOMModelDRM(void) : GModelData()
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
 * Constructs a COMPTEL DRM model from the information that is found in an
 * XML element. Please refer to the method GCOMModelDRM::read to learn more
 * about the information that is expected in the XML element.
 ***************************************************************************/
GCOMModelDRM::GCOMModelDRM(const GXmlElement& xml) : GModelData(xml)
{
    // Initialise members
    init_members();

    // Read XML
    read(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model COMPTEL DRM model.
 *
 * Constructs a COMPTEL DRM model by copying information from an existing
 * model. Note that the copy is a deep copy, so the original object can be
 * destroyed after the copy without any loss of information.
 ***************************************************************************/
GCOMModelDRM::GCOMModelDRM(const GCOMModelDRM& model) : GModelData(model)
{
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(model);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 *
 * Destroys a COMPTEL DRM model.
 ***************************************************************************/
GCOMModelDRM::~GCOMModelDRM(void)
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
 * @param[in] model COMPTEL DRM model.
 * @return COMPTEL DRM model
 *
 * Assigns the information from a COMPTEL DRM model to the actual object.
 * Note that a deep copy of the information is performed, so the original
 * object can be destroyed after the assignment without any loss of
 * information.
 ***************************************************************************/
GCOMModelDRM& GCOMModelDRM::operator=(const GCOMModelDRM& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelData::operator=(model);

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
 * @brief Clear instance
 *
 * Resets the object to a clean initial state. All information that resided
 * in the object will be lost.
 ***************************************************************************/
void GCOMModelDRM::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelData::free_members();
    this->GModel::free_members();

    // Initialise members
    this->GModel::init_members();
    this->GModelData::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 *
 * @return Pointer to deep copy of COMPTEL DRM model.
 *
 * Clone a COMPTEL DRM model. Cloning performs a deep copy of the
 * information, so the original object can be destroyed after cloning without
 * any loss of information.
 ***************************************************************************/
GCOMModelDRM* GCOMModelDRM::clone(void) const
{
    return new GCOMModelDRM(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 * @param[in] gradients Compute gradients?
 * @return Background model value.
 *
 * @exception GException::invalid_argument
 *            Observation is not a COMPTEL observation.
 *            Event is not a COMPTEL event bin.
 *
 * Evaluates the COMPTEL DRM model.
 ***************************************************************************/
double GCOMModelDRM::eval(const GEvent&       event,
                          const GObservation& obs,
                          const bool&         gradients) const
{
    // Extract COMPTEL observation
    const GCOMObservation* observation = dynamic_cast<const GCOMObservation*>(&obs);
    if (observation == NULL) {
        std::string cls = std::string(typeid(&obs).name());
        std::string msg = "Observation of type \""+cls+"\" is not a COMPTEL "
                          "observations. Please specify a COMPTEL observation "
                          "as argument.";
        throw GException::invalid_argument(G_EVAL, msg);
    }

    // Extract COMPTEL event bin
    const GCOMEventBin* bin = dynamic_cast<const GCOMEventBin*>(&event);
    if (bin == NULL) {
        std::string cls = std::string(typeid(&event).name());
        std::string msg = "Event of type \""+cls+"\" is  not a COMPTEL event. "
                          "Please specify a COMPTEL event as argument.";
        throw GException::invalid_argument(G_EVAL, msg);
    }

    // Initialise value
    double value = 0.0;

    // Optionally initialise gradients
    if (gradients) {
        m_norm.factor_gradient(0.0);
    }

    // Get bin index
    int index = bin->index();

    // Get bin size.
    double size = bin->size();
    
    // Continue only if bin size is positive
    if (size > 0.0) {

        // Get DRM model value
        value = m_drm.map().pixels()[index] / size;

        // Optionally compute gradients
        if (gradients) {

            // Compute partial derivative
            double grad = (m_norm.is_free()) ? value * m_norm.scale() : 0.0;

            // Set gradient
            m_norm.factor_gradient(grad);

        } // endif: gradient computation was requested

        // Compute model value
        value *= m_norm.value();

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
            std::cout << "*** ERROR: GCOMModelDRM::eval";
            std::cout << "(index=" << index << "):";
            std::cout << " NaN/Inf encountered";
            std::cout << " (value=" << value;
            std::cout << ")" << std::endl;
        }
        #endif

    } // endif: binsize was positive

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Return spatially integrated DRM model
 *
 * @param[in] obsEng Measured event energy.
 * @param[in] obsTime Measured event time.
 * @param[in] obs Observation.
 * @return Spatially integrated DRM model.
 *
 * @exception GException::feature_not_implemented
 *            Method is not implemented.
 *
 * Spatially integrates the DRM model for a given measured event energy and
 * event time.
 *
 * @todo Implement method.
 ***************************************************************************/
double GCOMModelDRM::npred(const GEnergy&      obsEng,
                           const GTime&        obsTime,
                           const GObservation& obs) const
{
    // Initialise result
    double npred = 0.0;

    // Method is not implemented
    std::string msg = "Spatial integration of DRM model is not implemented.";
    throw GException::feature_not_implemented(G_NRED, msg);

    // Return
    return npred;
}


/***********************************************************************//**
 * @brief Return simulated list of events
 *
 * @param[in] obs Observation.
 * @param[in] ran Random number generator.
 * @return COMPTEL event cube.
 *
 * @exception GException::feature_not_implemented
 *            Method is not implemented.
 *
 * Draws a sample of events from the COMPTEL DRM model using a Monte Carlo
 * simulation.
 *
 * @todo Implement method.
 ***************************************************************************/
GCOMEventCube* GCOMModelDRM::mc(const GObservation& obs, GRan& ran) const
{
    // Initialise new event cube
    GCOMEventCube* cube = new GCOMEventCube;

    // Method is not implemented
    std::string msg = "Monte Carlo simulation of DRM model is not implemented.";
    throw GException::feature_not_implemented(G_MC, msg);

    // Return
    return cube;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * Read the COMPTEL DRM model from an XML element.
 *
 * The model has the following XML file syntax:
 *
 *     <source name="Model" type="DRM" file="drm.fits" instrument="COM">
 *       <parameter name="Normalization" ../>
 *     </source>
 ***************************************************************************/
void GCOMModelDRM::read(const GXmlElement& xml)
{
    // Verify number of model parameters
    gammalib::xml_check_parnum(G_READ, xml, 1);

    // Get parameter pointers
    const GXmlElement* norm = gammalib::xml_get_par(G_READ, xml, m_norm.name());

    // Read parameters
    m_norm.read(*norm);

    // Read filename
    m_filename = gammalib::xml_file_expand(xml, xml.attribute("file"));

    // Load DRM
    m_drm.load(m_filename);

    // Read model attributes
    read_attributes(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element.
 *
 * Write the COMPTEL DRM model into an XML element.
 *
 * The model has the following XML file syntax:
 *
 *     <source name="Model" type="DRM" file="drm.fits" instrument="COM">
 *       <parameter name="Normalization" ../>
 *     </source>
 ***************************************************************************/
void GCOMModelDRM::write(GXmlElement& xml) const
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

    // If no source with corresponding name was found then append one.
    // Set also the type and the instrument.
    if (src == NULL) {
        src = xml.append("source");
        src->attribute("name", name());
        src->attribute("type", type());
        if (instruments().length() > 0) {
            src->attribute("instrument", instruments());
        }
    }

    // Verify model type
    gammalib::xml_check_type(G_WRITE, *src, type());

    // Get XML parameters
    GXmlElement* norm = gammalib::xml_need_par(G_WRITE, *src, m_norm.name());

    // Write parameters
    m_norm.write(*norm);

    // Set DRM file name
    src->attribute("file", gammalib::xml_file_reduce(*src, m_filename));

    // Write model attributes
    write_attributes(*src);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print model information
 *
 * @param[in] chatter Chattiness.
 * @return String containing model information.
 ***************************************************************************/
std::string GCOMModelDRM::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMModelDRM ===");

        // Append attributes
        result.append("\n"+print_attributes());

        // Append parameters
        result.append("\n"+gammalib::parformat("DRM file"));
        result.append(m_filename);
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
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCOMModelDRM::init_members(void)
{
    // Initialise Value
    m_norm.clear();
    m_norm.name("Normalization");
    m_norm.value(1.0);
    m_norm.scale(1.0);
    m_norm.range(0.0, 1.0e6);
    m_norm.gradient(0.0);
    m_norm.free();
    m_norm.has_grad(true);

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);

    // Initialise other members
    m_drm.clear();
    m_filename.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Model.
 ***************************************************************************/
void GCOMModelDRM::copy_members(const GCOMModelDRM& model)
{
    // Copy members
    m_norm     = model.m_norm;
    m_drm      = model.m_drm;
    m_filename = model.m_filename;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMModelDRM::free_members(void)
{
    // Return
    return;
}
