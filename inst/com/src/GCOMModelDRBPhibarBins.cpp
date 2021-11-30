/***************************************************************************
 * GCOMModelDRBPhibarBins.cpp - COMPTEL DRB model Phibar bin fitting class *
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
 * @file GCOMModelDRBPhibarBins.cpp
 * @brief COMPTEL DRB model Phibar bin fitting class implementation
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
#include "GEvent.hpp"
#include "GObservation.hpp"
#include "GXmlElement.hpp"
#include "GCOMModelDRBPhibarBins.hpp"
#include "GCOMObservation.hpp"
#include "GCOMEventBin.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GCOMModelDRBPhibarBins g_com_drb_bin_fitting_seed;
const GModelRegistry         g_com_drb_bin_fitting_registry(&g_com_drb_bin_fitting_seed);

/* __ Method name definitions ____________________________________________ */
#define G_EVAL  "GCOMModelDRBPhibarBins::eval(GEvent&, GObservation&, bool&)"
#define G_NRED             "GCOMModelDRBPhibarBins::npred(GEnergy&, GTime&, "\
                                                             "GObservation&)"
#define G_MC               "GCOMModelDRBPhibarBins::mc(GObservation&, GRan&)"
#define G_READ                   "GCOMModelDRBPhibarBins::read(GXmlElement&)"
#define G_WRITE                 "GCOMModelDRBPhibarBins::write(GXmlElement&)"

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
 * Constructs an empty COMPTEL DRB Phibar bin fitting model.
 ***************************************************************************/
GCOMModelDRBPhibarBins::GCOMModelDRBPhibarBins(void) : GModelData()
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
 * Constructs a COMPTEL DRB Phibar bin fitting model from the information
 * that is found in an XML element. Please refer to the method
 * GCOMModelDRBPhibarBins::read to learn more about the information that is
 * expected in the XML element.
 ***************************************************************************/
GCOMModelDRBPhibarBins::GCOMModelDRBPhibarBins(const GXmlElement& xml) :
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
 * @param[in] model COMPTEL DRB Phibar bin fitting model.
 *
 * Constructs a COMPTEL DRB Phibar bin fitting model by copying information
 * from an existing model. Note that the copy is a deep copy, so the original
 * object can be destroyed after the copy without any loss of information.
 ***************************************************************************/
GCOMModelDRBPhibarBins::GCOMModelDRBPhibarBins(const GCOMModelDRBPhibarBins& model) :
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
 * Destroys a COMPTEL DRB Phibar bin fitting model.
 ***************************************************************************/
GCOMModelDRBPhibarBins::~GCOMModelDRBPhibarBins(void)
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
 * @param[in] model COMPTEL DRB Phibar bin fitting model.
 * @return COMPTEL DRB Phibar bin fitting model
 *
 * Assigns the information from a COMPTEL DRB Phibar bin fitting model to the
 * actual object. Note that a deep copy of the information is performed, so
 * the original object can be destroyed after the assignment without any loss
 * of information.
 ***************************************************************************/
GCOMModelDRBPhibarBins& GCOMModelDRBPhibarBins::operator=(const GCOMModelDRBPhibarBins& model)
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
void GCOMModelDRBPhibarBins::clear(void)
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
 * @return Pointer to deep copy of COMPTEL DRB Phibar bin fitting model.
 *
 * Clone a COMPTEL DRB Phibar bin fitting model. Cloning performs a deep copy
 * of the information, so the original object can be destroyed after cloning
 * without any loss of information.
 ***************************************************************************/
GCOMModelDRBPhibarBins* GCOMModelDRBPhibarBins::clone(void) const
{
    return new GCOMModelDRBPhibarBins(*this);
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
 * @exception GException::invalid_value
 *            Model incompatible with data space.
 *
 * Evaluates the COMPTEL DRB Phibar bin fitting model.
 ***************************************************************************/
double GCOMModelDRBPhibarBins::eval(const GEvent&       event,
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

    // Get total number of Phibar layers and number of pixels
    int nphi = observation->drb().map().nmaps();
    int npix = observation->drb().map().npix();
    if (nphi != m_values.size()) {
        std::string msg = "Number of "+gammalib::str(nphi)+" Phibar layers "
                          "differs from number of "+gammalib::str(m_values.size())+
                          " DRB Phibar bin fitting model parameters. Please "
                          "specify a model that is compatible with the data "
                          "space.";
        throw GException::invalid_value(G_EVAL, msg);
    }
    if (npix == 0) {
        std::string msg = "There are no Chi/Psi pixels in the DRB model. "
                          "Please specify a model that contains Chi/Psi "
                          "pixels.";
        throw GException::invalid_value(G_EVAL, msg);
    }

    // Initialise value
    double value = 0.0;

    // Optionally initialise gradients
    if (gradients) {
        for (int i = 0; i < m_values.size(); ++i) {
            m_values[i].factor_gradient(0.0);
        }
    }

    // Get bin size.
    double size = bin->size();
    
    // Continue only if bin size is positive
    if (size > 0.0) {

        // Get bin index and Phibar layer
        int index = bin->index();
        int iphi  = index / npix;

        // Get DRB model value
        value = observation->drb().map().pixels()[index] / size;

        // Get scale factor
        double scale = m_values[iphi].value();

        // Optionally compute gradients
        if (gradients) {

            // Compute partial derivative
            double grad = (m_values[iphi].is_free())
                           ? value * m_values[iphi].scale() : 0.0;

            // Set gradient
            m_values[iphi].factor_gradient(grad);

        } // endif: gradient computation was requested

        // Compute background value
        value *= scale;

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
            std::cout << "*** ERROR: GCOMModelDRBPhibarBins::eval";
            std::cout << "(index=" << index << "):";
            std::cout << " NaN/Inf encountered";
            std::cout << " (iphi=" << iphi;
            std::cout << ", value=" << value;
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
double GCOMModelDRBPhibarBins::npred(const GEnergy&      obsEng,
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
 * Draws a sample of events from the COMPTEL DRB Phibar bin fitting model
 * using a Monte Carlo simulation.
 *
 * @todo Implement method.
 ***************************************************************************/
GCOMEventCube* GCOMModelDRBPhibarBins::mc(const GObservation& obs,
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
 * Read the COMPTEL DRB Phibar bin fitting model from an XML element.
 *
 * The model is composed of normalization values that define the scale
 * factors for all Phibar layers. The following XML file syntax is
 * expected:
 *
 *     <source name="Background" type="DRBPhibarBins" instrument="COM">
 *       <parameter name="Normalization" .../>
 *       ...
 *       <parameter name="Normalization" .../>
 *     </source>
 ***************************************************************************/
void GCOMModelDRBPhibarBins::read(const GXmlElement& xml)
{
    // Free space for bins
    m_values.clear();

    // Get number of bins from XML file
    int bins = xml.elements("parameter");

    // Loop over all bins
    for (int i = 0; i < bins; ++i) {

        // Allocate normalization parameter
        GModelPar normalization;
            
        // Read normalization parameter
        const GXmlElement* par = xml.element("parameter", i);
        normalization.read(*par);

        // Set parameter name
        std::string normalization_name = "Phibar normalization "+gammalib::str(i);

        // Set normalization attributes
        normalization.name(normalization_name);
        normalization.unit("");
        normalization.has_grad(true);

        // Push normalization parameter on list
        m_values.push_back(normalization);

    } // endfor: looped over bins

    // Read model attributes
    read_attributes(xml);

    // Set pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::invalid_value
 *            Existing XML element is not of required type
 *            Invalid number of model parameters or bins found in XML element.
 *
 * Write the COMPTEL DRB Phibar bin fitting model into an XML element.
 *
 * The model is composed of normalization values that define the scale
 * factors for all Phibar layers. The following XML file syntax is
 * expected:
 *
 *     <source name="Background" type="DRBPhibarBins" instrument="COM">
 *       <parameter name="Normalization" .../>
 *       ...
 *       <parameter name="Normalization" .../>
 *     </source>
 ***************************************************************************/
void GCOMModelDRBPhibarBins::write(GXmlElement& xml) const
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

    // Determine number of bins
    int bins = m_values.size();

    // If XML element has 0 parameters then append parameters
    if (src->elements() == 0) {
        for (int i = 0; i < bins; ++i) {
            src->append(GXmlElement("parameter"));
        }
    }

    // Verify that XML element has the required number of bins
    if (src->elements() != bins || src->elements("parameter") != bins) {
        std::string msg = "Invalid number of bins "+
                          gammalib::str(src->elements("parameter"))+
                          " found in XML file, but model requires exactly "+
                          gammalib::str(bins)+" bins.";
        throw GException::invalid_value(G_WRITE, msg);
    }

    // Loop over all bins
    for (int i = 0; i < bins; ++i) {

        // Get bin
        GXmlElement* par  = src->element("parameter", i);
        GModelPar    mpar = m_values[i];
        mpar.name("Normalization");
        mpar.write(*par);

    } // endfor: looped over bins

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
std::string GCOMModelDRBPhibarBins::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMModelDRBPhibarBins ===");

        // Append attributes
        result.append("\n"+print_attributes());

        // Append bins summary
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
void GCOMModelDRBPhibarBins::init_members(void)
{
    // Initialise members
    m_values.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Model.
 ***************************************************************************/
void GCOMModelDRBPhibarBins::copy_members(const GCOMModelDRBPhibarBins& model)
{
    // Copy members
    m_values = model.m_values;

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMModelDRBPhibarBins::free_members(void)
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
void GCOMModelDRBPhibarBins::set_pointers(void)
{
    // Clear parameter pointer(s)
    m_pars.clear();

    // Get number of bins
    int bins = m_values.size();

    // Set parameter pointers for all bins
    for (int i = 0; i < bins; ++i) {

        // Signal that values have gradients
        m_values[i].has_grad(true);

        // Set pointer
        m_pars.push_back(&(m_values[i]));

    } // endfor: looped over bins

    // Return
    return;
}
