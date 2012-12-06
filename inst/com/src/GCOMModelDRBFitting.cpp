/***************************************************************************
 *       GCOMModelDRBFitting.cpp  -  COMPTEL DRB model fitting class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelRegistry.hpp"
#include "GCOMModelDRBFitting.hpp"
#include "GCOMException.hpp"
#include "GCOMObservation.hpp"
#include "GCOMEventBin.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GCOMModelDRBFitting g_com_drb_fitting_seed;
const GModelRegistry      g_com_drb_fitting_registry(&g_com_drb_fitting_seed);

/* __ Method name definitions ____________________________________________ */
#define G_EVAL             "GCOMModelDRBFitting::eval(GEvent&,GObservation&)"
#define G_EVAL_GRADIENTS       "GCOMModelDRBFitting::eval_gradients(GEvent&,"\
                                                             "GObservation&)"
#define G_READ                      "GCOMModelDRBFitting::read(GXmlElement&)"
#define G_WRITE                    "GCOMModelDRBFitting::write(GXmlElement&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DUMP_MC                                  //!< Dump MC information


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
 * @brief Constructor
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
 * @return Background model value.
 *
 * @exception GCOMException::bad_observation_type
 *            Observation is not a COMPTEL observation.
 * @exception GCOMException::bad_event_type
 *            Event is not a COMPTEL event bin.
 *
 * Evaluates the COMPTEL DRB fitting model.
 ***************************************************************************/
double GCOMModelDRBFitting::eval(const GEvent&       event,
                                 const GObservation& obs) const
{
    // Extract COMPTEL observation
    const GCOMObservation* observation = dynamic_cast<const GCOMObservation*>(&obs);
    if (observation == NULL) {
        throw GCOMException::bad_observation_type(G_EVAL);
    }

    // Extract COMPTEL event bin
    const GCOMEventBin* bin = dynamic_cast<const GCOMEventBin*>(&event);
    if (bin == NULL) {
        throw GCOMException::bad_event_type(G_EVAL);
    }

    // Initialise value
    double value = 0.0;

    // Initialise scaling factor
    double scale = 1.0;

    // Get bin index
    int index = bin->index();

    // Get bin size.
    double size = bin->size();
    
    // Continue only if bin size is positive
    if (size > 0.0) {

        // Get DRB model value
        value = observation->drb().pixels()[index] / size;

        // If model is a scaling factor then use the single parameter as such
        if (m_scale) {
            scale = m_values[0].real_value();
        }

        // ... otherwise perform a linear interpolation
        else {

            // Get Phibar value
            double phibar = bin->dir().phibar();

            // Update evaluation cache
            update_cache();

            // Set node array for linear interpolation
            m_nodes.set_value(phibar);

            // Get scale factor
            scale = m_values[m_nodes.inx_left()].real_value()  * m_nodes.wgt_left() +
                    m_values[m_nodes.inx_right()].real_value() * m_nodes.wgt_right();

        } // endelse: performed linear interpolation

        // Compute background value
        value *= scale;

    } // endif: binsize was positive

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (isnotanumber(value) || isinfinite(value)) {
        std::cout << "*** ERROR: GCOMModelDRBFitting::eval";
        std::cout << "(index=" << index << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", scale=" << scale;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 * @return Background model value.
 *
 * @exception GCOMException::bad_observation_type
 *            Observation is not a COMPTEL observation.
 * @exception GCOMException::bad_event_type
 *            Event is not a COMPTEL event bin.
 *
 * Evaluates the COMPTEL DRB fitting model value and sets the parameter
 * gradients.
 ***************************************************************************/
double GCOMModelDRBFitting::eval_gradients(const GEvent&       event,
                                           const GObservation& obs) const
{
    // Extract COMPTEL observation
    const GCOMObservation* observation = dynamic_cast<const GCOMObservation*>(&obs);
    if (observation == NULL) {
        throw GCOMException::bad_observation_type(G_EVAL_GRADIENTS);
    }

    // Extract COMPTEL event bin
    const GCOMEventBin* bin = dynamic_cast<const GCOMEventBin*>(&event);
    if (bin == NULL) {
        throw GCOMException::bad_event_type(G_EVAL_GRADIENTS);
    }

    // Initialise value and gradients
    double value = 0.0;
    for (int i = 0; i < m_values.size(); ++i) {
        const_cast<GCOMModelDRBFitting*>(this)->m_values[i].gradient(0.0);
    }

    // Initialise scaling factor
    double scale = 1.0;

    // Get bin index
    int index = bin->index();

    // Get bin size.
    double size = bin->size();
    
    // Continue only if bin size is positive
    if (size > 0.0) {

        // Get DRB model value
        value = observation->drb().pixels()[index] / size;


        // If model is a scaling factor then use the single parameter as such
        if (m_scale) {

            // Get scale factor
            scale = m_values[0].real_value();
        
            // Compute partial derivative
            double grad = (m_values[0].isfree()) ? value * m_values[0].scale() : 0.0;

            // Set gradient (circumvent const correctness)
            const_cast<GCOMModelDRBFitting*>(this)->m_values[0].gradient(grad);

        }

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
            scale = m_values[inx_left].real_value()  * wgt_left +
                    m_values[inx_right].real_value() * wgt_right;

            // Gradient for left node
            if (m_values[inx_left].isfree()) {
                double grad = wgt_left * value * m_values[inx_left].scale();
                const_cast<GCOMModelDRBFitting*>(this)->m_values[inx_left].gradient(grad);
            }

            // Gradient for right node
            if (m_values[inx_right].isfree()) {
                double grad = wgt_right * value * m_values[inx_right].scale();
                const_cast<GCOMModelDRBFitting*>(this)->m_values[inx_right].gradient(grad);
            }

        } // endelse: performed linear interpolation

        // Compute background value
        value *= scale;

    } // endif: binsize was positive

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (isnotanumber(value) || isinfinite(value)) {
        std::cout << "*** ERROR: GCOMModelDRBFitting::eval_gradients";
        std::cout << "(index=" << index << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", scale=" << scale;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Return spatially integrated data model
 *
 * @param[in] obsEng Measured event energy.
 * @param[in] obsTime Measured event time.
 * @param[in] obs Observation.
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

    // Return
    return npred;
}


/***********************************************************************//**
 * @brief Return simulated list of events
 *
 * @param[in] obs Observation.
 * @param[in] ran Random number generator.
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

    // Return
    return cube;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter name found in XML element.
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
        std::string message = "DRB fitting model requires at least one"
                              " Phibar node.";
        throw GException::model_invalid_parnum(G_READ, xml, message);
    }

    // Loop over all nodes
    for (int i = 0; i < nodes; ++i) {

        // Allocate node parameters
        GModelPar phibar;
        GModelPar normalization;
            
        // Get node
        GXmlElement* node = static_cast<GXmlElement*>(xml.element("node", i));

        // Verify that node XML element has exactly 2 parameters
        if (node->elements() != 2 || node->elements("parameter") != 2) {
            throw GException::model_invalid_parnum(G_READ, xml,
                  "Node requires exactly 2 parameters.");
        }

        // Extract node parameters
        int npar[] = {0, 0};
        for (int k = 0; k < 2; ++k) {

            // Get parameter element
            GXmlElement* par = static_cast<GXmlElement*>(node->element("parameter", k));

            // Handle energy
            if (par->attribute("name") == "Phibar") {
                phibar.read(*par);
                npar[0]++;
            }

            // Handle intensity
            else if (par->attribute("name") == "Normalization") {
                normalization.read(*par);
                npar[1]++;
            }

        } // endfor: looped over parameters

        // Verify that all parameters were found
        if (npar[0] != 1 || npar[1] != 1) {
            throw GException::model_invalid_parnames(G_READ, xml,
                  "Require \"Phibar\" and \"Normalization\" parameters.");
        }

        // Set parameter names
        std::string phibar_name        = "Phibar layer "+str(i);
        std::string normalization_name = "Scale factor "+str(i);

        // Set Phibar attributes
        phibar.name(phibar_name);
        phibar.unit("deg");
        phibar.hasgrad(false);

        // Set normalization attributes
        normalization.name(normalization_name);
        normalization.unit("");
        normalization.hasgrad(true);

        // Push node parameters on list
        m_phibars.push_back(phibar);
        m_values.push_back(normalization);

    } // endfor: looped over nodes

    // Set model name
    name(xml.attribute("name"));

    // Set instruments
    instruments(xml.attribute("instrument"));

    // Set observation identifiers
    ids(xml.attribute("id"));

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
        GXmlElement* element = static_cast<GXmlElement*>(xml.element("source", k));
        if (element->attribute("name") == name()) {
            src = element;
            break;
        }
    }

    // If no source with corresponding name was found then append one.
    // Set also the type and the instrument.
    if (src == NULL) {
        src = new GXmlElement("source");
        src->attribute("name", name());
        src->attribute("type", type());
        if (instruments().length() > 0) {
            src->attribute("instrument", instruments());
        }
        xml.append(src);
    }

    // Verify model type
    if (src->attribute("type") != type()) {
        throw GException::model_invalid_spectral(G_WRITE, src->attribute("type"),
              "DRB fitting model is not of type \""+type()+"\".");
    }

    // Determine number of nodes
    int nodes = m_phibars.size();

    // If XML element has 0 nodes then append nodes
    if (src->elements() == 0) {
        for (int i = 0; i < nodes; ++i) {
            src->append(new GXmlElement("node"));
        }
    }

    // Verify that XML element has the required number of nodes
    if (src->elements() != nodes || src->elements("node") != nodes) {
        std::string message = "DRB fitting model requires exactly " +
                              str(nodes) + " nodes.";
        throw GException::model_invalid_parnum(G_WRITE, *src, message);
    }

    // Loop over all nodes
    for (int i = 0; i < nodes; ++i) {

        // Get node
        GXmlElement* node = static_cast<GXmlElement*>(src->element("node", i));

        // If XML element has 0 leafs then append energy and intensity
        // element
        if (node->elements() == 0) {
            node->append(new GXmlElement("parameter name=\"Phibar\""));
            node->append(new GXmlElement("parameter name=\"Normalization\""));
        }

        // Verify that node XML element has exactly 2 parameters
        if (node->elements() != 2 || node->elements("parameter") != 2) {
            throw GException::model_invalid_parnum(G_WRITE, *src,
                  "Node requires exactly 2 parameters.");
        }

        // Set or update model parameter attributes
        int npar[] = {0, 0};
        for (int k = 0; k < 2; ++k) {

            // Get parameter element
            GXmlElement* par = static_cast<GXmlElement*>(node->element("parameter", k));

            // Handle prefactor
            if (par->attribute("name") == "Phibar") {
                npar[0]++;
                m_phibars[i].write(*par);
            }

            // Handle index
            else if (par->attribute("name") == "Normalization") {
                npar[1]++;
                m_values[i].write(*par);
            }

        } // endfor: looped over parameters

        // Check of all required parameters are present
        if (npar[0] != 1 || npar[1] != 1) {
            throw GException::model_invalid_parnames(G_WRITE, *src,
                  "Require \"Phibar\" and \"Normalization\" parameters.");
        }

    } // endfor: looped over nodes

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print model information
 *
 * @todo Document method.
 ***************************************************************************/
std::string GCOMModelDRBFitting::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GCOMModelDRBFitting ===");

    // Append attributes
    result.append("\n"+print_attributes());

    // Append node summary
    result.append("\n"+parformat("Number of nodes")+str(m_phibars.size()));
    result.append("\n"+parformat("Number of parameters")+str(size()));
    for (int i = 0; i < size(); ++i) {
        result.append("\n"+m_pars[i]->print());
    }

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
        m_phibars[i].hasgrad(false);

        // Signal that values have gradients
        m_values[i].hasgrad(true);

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
            double phibar = m_phibars[i].real_value();
            m_nodes.append(phibar);
            m_old_phibars.push_back(phibar);
            if (m_phibars[i].isfree()) {
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
            double phibar = m_phibars[i].real_value();
            if (phibar != m_old_phibars[i]) {
                m_nodes[i]       = phibar;
                m_old_phibars[i] = phibar;
            }
        }

    } // endif: there were at least two nodes

    // Return
    return;
}
