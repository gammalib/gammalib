/***************************************************************************
 *        GCOMModelDRBFitting.cpp - COMPTEL DRB model fitting class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2021 by Juergen Knoedlseder                         *
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
 * @file GCOMModelDRBFitting.cpp
 * @brief COMPTEL DRB model fitting class implementation
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
#include "GCOMModelDRBFitting.hpp"
#include "GCOMObservation.hpp"
#include "GCOMEventBin.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GCOMModelDRBFitting g_com_drb_fitting_seed;
const GModelRegistry      g_com_drb_fitting_registry(&g_com_drb_fitting_seed);

/* __ Method name definitions ____________________________________________ */
#define G_EVAL     "GCOMModelDRBFitting::eval(GEvent&, GObservation&, bool&)"
#define G_NRED  "GCOMModelDRBFitting::npred(GEnergy&, GTime&, GObservation&)"
#define G_MC                  "GCOMModelDRBFitting::mc(GObservation&, GRan&)"
#define G_READ                      "GCOMModelDRBFitting::read(GXmlElement&)"
#define G_WRITE                    "GCOMModelDRBFitting::write(GXmlElement&)"

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
 * Constructs an empty COMPTEL DRB fitting model.
 ***************************************************************************/
GCOMModelDRBFitting::GCOMModelDRBFitting(void) : GModelData()
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
 * Constructs a COMPTEL DRB fitting model from the information that is found
 * in an XML element. Please refer to the method GCOMModelDRBFitting::read
 * to learn more about the information that is expected in the XML element.
 ***************************************************************************/
GCOMModelDRBFitting::GCOMModelDRBFitting(const GXmlElement& xml) :
                     GModelData(xml)
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
 * @param[in] model COMPTEL DRB fitting model.
 *
 * Constructs a COMPTEL DRB fitting model by copying information from an
 * existing model. Note that the copy is a deep copy, so the original object
 * can be destroyed after the copy without any loss of information.
 ***************************************************************************/
GCOMModelDRBFitting::GCOMModelDRBFitting(const GCOMModelDRBFitting& model) :
                     GModelData(model)
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
 * Destroys a COMPTEL DRB fitting model.
 ***************************************************************************/
GCOMModelDRBFitting::~GCOMModelDRBFitting(void)
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
 * @param[in] model COMPTEL DRB fitting model.
 * @return COMPTEL DRB fitting model
 *
 * Assigns the information from a COMPTEL DRB fitting model to the actual
 * object. Note that a deep copy of the information is performed, so the
 * original object can be destroyed after the assignment without any loss of
 * information.
 ***************************************************************************/
GCOMModelDRBFitting& GCOMModelDRBFitting::operator=(const GCOMModelDRBFitting& model)
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
void GCOMModelDRBFitting::clear(void)
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
 * @return Pointer to deep copy of COMPTEL DRB fitting model.
 *
 * Clone a COMPTEL DRB fitting model. Cloning performs a deep copy of the
 * information, so the original object can be destroyed after cloning without
 * any loss of information.
 ***************************************************************************/
GCOMModelDRBFitting* GCOMModelDRBFitting::clone(void) const
{
    return new GCOMModelDRBFitting(*this);
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
 * Evaluates the COMPTEL DRB fitting model.
 ***************************************************************************/
double GCOMModelDRBFitting::eval(const GEvent&       event,
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
        for (int i = 0; i < m_values.size(); ++i) {
            m_values[i].factor_gradient(0.0);
        }
    }

    // Get bin index
    int index = bin->index();

    // Get bin size.
    double size = bin->size();
    
    // Continue only if bin size is positive
    if (size > 0.0) {

        // Initialise scaling factor
        double scale = 1.0;

        // Get DRB model value
        value = observation->drb().map().pixels()[index] / size;

        // If the model has a single scaling factor then use the first
        // parameter as scaling factor of the model
        if (m_scale) {

            // Get scale factor
            scale = m_values[0].value();

            // Optionally compute gradients
            if (gradients) {

                // Compute partial derivative
                double grad = (m_values[0].is_free())
                              ? value * m_values[0].scale() : 0.0;

                // Set gradient
                m_values[0].factor_gradient(grad);

            } // endif: gradient computation was requested

        } // endif: model is a scaling factor

        // ... otherwise perform a linear interpolation
        else {

            // Get Phibar value
            double phibar = bin->dir().phibar();

            // Update evaluation cache
            update_cache();

            // Set node array for linear interpolation
            m_nodes.set_value(phibar);

            // Get indices and weights for interpolation
            int    inx_left  = m_nodes.inx_left();
            int    inx_right = m_nodes.inx_right();
            double wgt_left  = m_nodes.wgt_left();
            double wgt_right = m_nodes.wgt_right();

            // Get scale factor
            scale = m_values[inx_left].value()  * wgt_left +
                    m_values[inx_right].value() * wgt_right;

            // Optionally compute gradients
            if (gradients) {

                // Compute partial derivatives
                double g_left  = (m_values[inx_left].is_free())
                                 ? wgt_left * value * m_values[inx_left].scale()
                                 : 0.0;
                double g_right = (m_values[inx_right].is_free())
                                 ? wgt_right * value * m_values[inx_right].scale()
                                 : 0.0;

                // Set gradients
                m_values[inx_left].factor_gradient(g_left);
                m_values[inx_right].factor_gradient(g_right);

            } // endif: gradient computation was requested

        } // endelse: performed linear interpolation

        // Compute background value
        value *= scale;

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
            std::cout << "*** ERROR: GCOMModelDRBFitting::eval";
            std::cout << "(index=" << index << "):";
            std::cout << " NaN/Inf encountered";
            std::cout << " (value=" << value;
            std::cout << ", scale=" << scale;
            std::cout << ")" << std::endl;
        }
        #endif

    } // endif: binsize was positive

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Return spatially integrated data model
 *
 * @param[in] obsEng Measured event energy.
 * @param[in] obsTime Measured event time.
 * @param[in] obs Observation.
 * @return Spatially integrated data model.
 *
 * @exception GException::feature_not_implemented
 *            Method is not implemented.
 *
 * Spatially integrates the data model for a given measured event energy and
 * event time.
 *
 * @todo Implement method.
 ***************************************************************************/
double GCOMModelDRBFitting::npred(const GEnergy&      obsEng,
                                  const GTime&        obsTime,
                                  const GObservation& obs) const
{
    // Initialise result
    double npred = 0.0;

    // Method is not implemented
    std::string msg = "Spatial integration of data model is not implemented.";
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
 * Draws a sample of events from the COMPTEL DRB fitting model using a Monte
 * Carlo simulation.
 *
 * @todo Implement method.
 ***************************************************************************/
GCOMEventCube* GCOMModelDRBFitting::mc(const GObservation& obs, 
                                       GRan&               ran) const
{
    // Initialise new event cube
    GCOMEventCube* cube = new GCOMEventCube;

    // Method is not implemented
    std::string msg = "Monte Carlo simulation of data model is not implemented.";
    throw GException::feature_not_implemented(G_MC, msg);

    // Return
    return cube;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::invalid_value
 *            Model definition requires at least one node.
 *
 * Read the COMPTEL DRB fitting model from an XML element.
 *
 * The model is composed of nodes that define the normalization value as
 * function of Phibar value. The following XML file syntax is expected:
 *
 *     <source name="Background" type="DRBFitting" instrument="COM">
 *       <node>
 *         <parameter name="Phibar"        .../>
 *         <parameter name="Normalization" .../>
 *       </node>
 *       ...
 *       <node>
 *         <parameter name="Phibar"        .../>
 *         <parameter name="Normalization" .../>
 *       </node>
 *     </source>
 ***************************************************************************/
void GCOMModelDRBFitting::read(const GXmlElement& xml)
{
    // Free space for nodes
    m_phibars.clear();
    m_values.clear();

    // Get number of nodes from XML file
    int nodes = xml.elements("node");

    // Throw an error if there are no nodes
    if (nodes < 1) {
        std::string msg = "DRB fitting model requires at least one Phibar "
                          "node. Please correct XML format.";
        throw GException::invalid_value(G_READ, msg);
    }

    // Loop over all nodes
    for (int i = 0; i < nodes; ++i) {

        // Allocate node parameters
        GModelPar phibar;
        GModelPar normalization;
            
        // Get node
        const GXmlElement* node = xml.element("node", i);

        // Read Phibar parameter
        const GXmlElement* par = gammalib::xml_get_par(G_READ, *node, "Phibar");
        phibar.read(*par);

        // Read Normalization parameter
        par = gammalib::xml_get_par(G_READ, *node, "Normalization");
        normalization.read(*par);

        // Set parameter names
        std::string phibar_name        = "Phibar layer "+gammalib::str(i);
        std::string normalization_name = "Scale factor "+gammalib::str(i);

        // Set Phibar attributes
        phibar.name(phibar_name);
        phibar.unit("deg");
        phibar.has_grad(false);

        // Set normalization attributes
        normalization.name(normalization_name);
        normalization.unit("");
        normalization.has_grad(true);

        // Push node parameters on list
        m_phibars.push_back(phibar);
        m_values.push_back(normalization);

    } // endfor: looped over nodes

    // Read model attributes
    read_attributes(xml);

    // Set pointers
    set_pointers();

    // Set evluation cache
    set_cache();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::model_invalid_spectral
 *            Existing XML element is not of required type
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters or nodes found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Write the COMPTEL DRB fitting model into an XML element.
 *
 * The model is composed of nodes that define the normalization value as
 * function of Phibar value. The following XML file syntax is expected:
 *
 *     <source name="Background" type="DRBFitting" instrument="COM">
 *       <node>
 *         <parameter name="Phibar"        .../>
 *         <parameter name="Normalization" .../>
 *       </node>
 *       ...
 *       <node>
 *         <parameter name="Phibar"        .../>
 *         <parameter name="Normalization" .../>
 *       </node>
 *     </source>
 ***************************************************************************/
void GCOMModelDRBFitting::write(GXmlElement& xml) const
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
    if (src->attribute("type") != type()) {
        std::string msg = "Invalid model type \""+src->attribute("type")+"\" "
                          "found in XML file. Model type \""+type()+"\" "
                          "expected.";
        throw GException::invalid_value(G_WRITE, msg);
    }

    // Determine number of nodes
    int nodes = m_phibars.size();

    // If XML element has 0 nodes then append nodes
    if (src->elements() == 0) {
        for (int i = 0; i < nodes; ++i) {
            src->append(GXmlElement("node"));
        }
    }

    // Verify that XML element has the required number of nodes
    if (src->elements() != nodes || src->elements("node") != nodes) {
        std::string msg = "Invalid number of nodes "+
                          gammalib::str(src->elements("node"))+
                          " found in XML file, but model requires exactly "+
                          gammalib::str(nodes)+" nodes.";
        throw GException::invalid_value(G_WRITE, msg);
    }

    // Loop over all nodes
    for (int i = 0; i < nodes; ++i) {

        // Get node
        GXmlElement* node = src->element("node", i);

        // Write Phibar parameter
        GXmlElement* par  = gammalib::xml_need_par(G_WRITE, *node, "Phibar");
        GModelPar    mpar = m_phibars[i];
        mpar.name("Phibar");
        mpar.write(*par);

        // Write Normalization parameter
        par  = gammalib::xml_need_par(G_WRITE, *node, "Normalization");
        mpar = m_values[i];
        mpar.name("Normalization");
        mpar.write(*par);

    } // endfor: looped over nodes

    // Write model attributes
    write_attributes(*src);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print model information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing model information.
 ***************************************************************************/
std::string GCOMModelDRBFitting::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMModelDRBFitting ===");

        // Append attributes
        result.append("\n"+print_attributes());

        // Append node summary
        result.append("\n"+gammalib::parformat("Number of nodes"));
        result.append(gammalib::str(m_phibars.size()));
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
void GCOMModelDRBFitting::init_members(void)
{
    // Initialise members
    m_phibars.clear();
    m_values.clear();

    // Initialise cache
    m_scale = true;
    m_fixed = true;
    m_old_phibars.clear();
    m_nodes.clear();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Model.
 ***************************************************************************/
void GCOMModelDRBFitting::copy_members(const GCOMModelDRBFitting& model)
{
    // Copy members
    m_phibars = model.m_phibars;
    m_values  = model.m_values;

    // Copy cache
    m_scale       = model.m_scale;
    m_fixed       = model.m_fixed;
    m_old_phibars = model.m_old_phibars;
    m_nodes       = model.m_nodes;

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMModelDRBFitting::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set pointers
 *
 * Set pointers to all model parameters. The pointers are stored in a vector
 * that is member of the GModel base class.
 ***************************************************************************/
void GCOMModelDRBFitting::set_pointers(void)
{
    // Clear parameter pointer(s)
    m_pars.clear();

    // Get number of nodes
    int nodes = m_phibars.size();

    // Set parameter pointers for all nodes
    for (int i = 0; i < nodes; ++i) {

        // Signal that Phibar values have no gradients
        m_phibars[i].has_grad(false);

        // Signal that values have gradients
        m_values[i].has_grad(true);

        // Set pointer
        m_pars.push_back(&(m_phibars[i]));
        m_pars.push_back(&(m_values[i]));

    } // endfor: looped over nodes

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set evaluation cache
 *
 * Set the evaluation cache for fast computations. The evaluation cache
 * consists of the Phibar values stored in a node array. The evaluation
 * cache is only filled if at least 2 nodes are available. If only a single
 * node is available, the model is considered as a simple scaling factor.
 ***************************************************************************/
void GCOMModelDRBFitting::set_cache(void) const
{
    // Clear any existing values
    m_old_phibars.clear();
    m_nodes.clear();
    m_fixed = true;
    m_scale = true; // Signals that model is a scale factor

    // Determine number of parameters
    int num = m_phibars.size();

    // Continue only if we have at least two nodes
    if (num >= 2) {

        // Signal that model is not a scale factor
        m_scale = false;

        // Set Phibar node array. Signal if one of the Phibar values is free
        for (int i = 0; i < num; ++i) {
            double phibar = m_phibars[i].value();
            m_nodes.append(phibar);
            m_old_phibars.push_back(phibar);
            if (m_phibars[i].is_free()) {
                m_fixed = false;
            }
        }

    } // endif: there were at least two nodes

    // Return
    return;
}


/***********************************************************************//**
 * @brief Update evaluation cache
 *
 * Updates the evaluation cache by computing the 
 * changed.
 ***************************************************************************/
void GCOMModelDRBFitting::update_cache(void) const
{
    // Determine number of parameters
    int num = m_phibars.size();

    // Continue only if we have at least two nodes
    if (num >= 2) {
    
        // Update Phibar values
        for (int i = 0; i < num; ++i) {
            double phibar = m_phibars[i].value();
            if (phibar != m_old_phibars[i]) {
                m_nodes[i]       = phibar;
                m_old_phibars[i] = phibar;
            }
        }

    } // endif: there were at least two nodes

    // Return
    return;
}
