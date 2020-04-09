/***************************************************************************
 *          GSPIModelDataSpace.cpp - INTEGRAL/SPI data space model         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GSPIModelDataSpace.cpp
 * @brief INTEGRAL/SPI data space model implementation
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
#include "GModelPar.hpp"
#include "GSPIModelDataSpace.hpp"
#include "GSPIObservation.hpp"
#include "GSPIEventBin.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GSPIModelDataSpace g_spi_data_space_seed;
const GModelRegistry     g_spi_data_space_registry(&g_spi_data_space_seed);

/* __ Method name definitions ____________________________________________ */
#define G_CONSTRUCTOR               "GSPIModelDataSpace::GSPIModelDataSpace("\
                        "GSPIObservation&, std::string&, std::string&, int&)"
#define G_EVAL      "GSPIModelDataSpace::eval(GEvent&, GObservation&, bool&)"
#define G_READ                       "GSPIModelDataSpace::read(GXmlElement&)"
#define G_WRITE                     "GSPIModelDataSpace::write(GXmlElement&)"
#define G_SETUP_MODEL         "GSPIModelDataSpace::setup_model(Observation&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DUMP_MC                                  //!< Dump MC information
//#define G_DEBUG_SETUP_MODEL                   //!< Debug setup model method
//#define G_DEBUG_SETUP_MODEL_MAP           //!< Debug setup of parameter map


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Constructs an empty INTEGRAL/SPI data space model.
 ***************************************************************************/
GSPIModelDataSpace::GSPIModelDataSpace(void) : GModelData()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Model and method constructor
 *
 * @param[in] obs INTEGRAL/SPI observation.
 * @param[in] name Model name.
 * @param[in] method Fit method.
 * @param[in] index Model index.
 *
 * @exception GException::out_of_range
 *            Model index out of range.
 *
 * Constructs an INTEGRAL/SPI data space model by setting the fit method
 * and the model index.
 *
 * The constructor uses the INTEGRAL/SPI observation to setup the model
 * parameters and the mapping of parameters to the data space.
 ***************************************************************************/
GSPIModelDataSpace::GSPIModelDataSpace(const GSPIObservation& obs,
                                       const std::string&     name,
                                       const std::string&     method,
                                       const int&             index) :
                    GModelData()
{
    // Initialise members
    init_members();

    // Set members
    m_name   = name;
    m_method = method;
    m_index  = index;

    // Setup model
    setup_model(obs);

    // Check whether the index is valid. We can only do this if the
    // observation contains an event cube
    GSPIEventCube* cube = dynamic_cast<GSPIEventCube*>(m_obs->events());
    if (cube != NULL) {
        int num_models = cube->models();
        if (index < 0 || index >= num_models) {
            throw GException::out_of_range(G_CONSTRUCTOR, "Invalid model index",
                                           index, num_models);
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs an INTEGRAL/SPI data space model from the information that is
 * found in an XML element. Please refer to the read() method to learn more
 * about the information that is expected in the XML element.
 ***************************************************************************/
GSPIModelDataSpace::GSPIModelDataSpace(const GXmlElement& xml) :
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
 * @param[in] model INTEGRAL/SPI data space model.
 *
 * Constructs an INTEGRAL/SPI data space model by copying information from an
 * existing model. Note that the copy is a deep copy, so the original object
 * can be destroyed after the copy without any loss of information.
 ***************************************************************************/
GSPIModelDataSpace::GSPIModelDataSpace(const GSPIModelDataSpace& model) :
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
 * Destroys an INTEGRAL/SPI data space model.
 ***************************************************************************/
GSPIModelDataSpace::~GSPIModelDataSpace(void)
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
 * @param[in] model INTEGRAL/SPI data space model.
 * @return INTEGRAL/SPI data space model
 *
 * Assigns the information from an INTEGRAL/SPI data space model to the
 * class instance. Note that a deep copy of the information is performed, so
 * the @p model instance can be destroyed after the assignment without any
 * loss of information.
 ***************************************************************************/
GSPIModelDataSpace& GSPIModelDataSpace::operator=(const GSPIModelDataSpace& model)
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
void GSPIModelDataSpace::clear(void)
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
 * @return Pointer to deep copy of INTEGRAL/SPI data space model.
 *
 * Clone an INTEGRAL/SPI data space model. Cloning performs a deep copy of
 * the information, so the original object can be destroyed after cloning
 * without any loss of information.
 ***************************************************************************/
GSPIModelDataSpace* GSPIModelDataSpace::clone(void) const
{
    return new GSPIModelDataSpace(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 * @param[in] gradients Compute gradients?
 * @return Model value.
 *
 * @exception GException::invalid_argument
 *            Observation is not an INTEGRAL/SPI observation.
 *            Event is not an INTEGRAL/SPI event bin.
 *
 * Evaluates the value of the INTEGRAL/SPI data space model for a given
 * event and a given observation.
 ***************************************************************************/
double GSPIModelDataSpace::eval(const GEvent&       event,
                                const GObservation& obs,
                                const bool&         gradients) const
{
    // Reset parameter indices
    m_has_eval_inx = false;
    m_eval_inx.clear();

    // Extract INTEGRAL/SPI event bin
    const GSPIEventBin* bin = dynamic_cast<const GSPIEventBin*>(&event);
    if (bin == NULL) {
        std::string cls = std::string(typeid(&event).name());
        std::string msg = "Event of type \""+cls+"\" is  not an INTEGRAL/SPI "
                          "event bin. Please specify an INTEGRAL/SPI event "
                          "bin as argument.";
        throw GException::invalid_argument(G_EVAL, msg);
    }

    // Setup model
    setup_model(obs);

    // Initialise value
    double value = 0.0;

    // Continue only if model was initialised
    if (m_index >= 0 && m_parameters.size() > 0) {

        // Get bin size
        double size = bin->size();

        // Continue only if bin size is positive
        if (size > 0.0) {

            // Get bin index in data space
            int index = bin->index();

            // Get parameter index
            int ipar = m_map[index];

            // Get relevant model parameter
            const GModelPar* par = &(m_parameters[ipar]);

            // Get model value divided by size (note that the size will be
            // multiplied-in later, hence here we need the model value divided
            // by size)
            value = bin->model(m_index) / size;

            // Get model scale factor
            double scale = par->value();

            // Optionally compute gradients
            if (gradients) {

                // Compute partial derivative
                double grad = (par->is_free()) ? value * par->scale() : 0.0;

                // Set gradient
                par->factor_gradient(grad);

                // Signal that gradient for this parameter was set
                m_has_eval_inx = true;
                m_eval_inx.push_back(ipar);

            } // endif: gradient computation was requested

            // Compute scaled model value
            value *= scale;

            // Compile option: Check for NaN/Inf
            #if defined(G_NAN_CHECK)
            if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
                std::cout << "*** ERROR: GSPIModelDataSpace::eval";
                std::cout << "(index=" << index << "):";
                std::cout << " NaN/Inf encountered";
                std::cout << " (value=" << value;
                std::cout << ", scale=" << scale;
                std::cout << ")" << std::endl;
            }
            #endif

        } // endif: binsize was positive

    } // endif: model was initialised

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
double GSPIModelDataSpace::npred(const GEnergy&      obsEng,
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
 * Draws a sample of events from the INTEGRAL/SPI data space model using a
 * Monte Carlo simulation.
 *
 * @todo Implement method.
 ***************************************************************************/
GSPIEventCube* GSPIModelDataSpace::mc(const GObservation& obs,
                                      GRan&               ran) const
{
    // Initialise new event cube
    GSPIEventCube* cube = new GSPIEventCube;

    // Return
    return cube;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * Read the INTEGRAL/SPI data space model from an XML element.
 *
 * The model is composed of a list of scale parameters. The following XML
 * file format is expected:
 *
 *     <source name="Background" type="DataSpace" method="orbit,dete" instrument="SPI">
 *       <parameter name="REV0019_DETID00"        .../>
 *       <parameter name="REV0019_DETID01"        .../>
 *       ...
 *     </source>
 ***************************************************************************/
void GSPIModelDataSpace::read(const GXmlElement& xml)
{
    // Clear instance
    clear();

    // Set model method
    m_method = xml.attribute("method");

    // Get number of parameters from XML file
    int npars = xml.elements("parameter");

    // Loop over all parameters
    for (int i = 0; i < npars; ++i) {

        // Allocate parameter
        GModelPar parameter;

        // Get XML element
        const GXmlElement* par = xml.element("parameter", i);

        // Read parameter from XML element
        parameter.read(*par);

        // Push parameter on list
        m_parameters.push_back(parameter);

    } // endfor: looped over parameters

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
 *            Invalid number of model parameters or nodes found in XML element.
 *
 * Write the INTEGRAL/SPI data space model into an XML element.
 *
 * The model is composed of a list of scale parameters. The following XML
 * file structure will be written:
 *
 *     <source name="Background" type="DataSpace" method="orbit,dete" instrument="SPI">
 *       <parameter name="REV0019_DETID00"        .../>
 *       <parameter name="REV0019_DETID01"        .../>
 *       ...
 *     </source>
 ***************************************************************************/
void GSPIModelDataSpace::write(GXmlElement& xml) const
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

    // Determine number of parameters
    int npars = m_parameters.size();

    // If XML element has 0 parameters then append parameters
    if (src->elements() == 0) {
        for (int i = 0; i < npars; ++i) {
            src->append(GXmlElement("parameter"));
        }
    }

    // Verify that XML element has the required number of nodes
    if (src->elements() != npars || src->elements("parameter") != npars) {
        std::string msg = "Invalid number of parameters "+
                          gammalib::str(src->elements("parameter"))+
                          " found in XML file, but model requires exactly "+
                          gammalib::str(npars)+" parameters.";
        throw GException::invalid_value(G_WRITE, msg);
    }

    // Loop over all parameters
    for (int i = 0; i < npars; ++i) {

        // Get XML element
        GXmlElement* par = src->element("parameter", i);

        // Write parameter into XML element
        m_parameters[i].write(*par);

    } // endfor: looped over parameters

    // Write model attributes
    write_attributes(*src);

    // Write method
    xml.attribute("method", m_method);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print model information
 *
 * @param[in] chatter Chattiness.
 * @return String containing model information.
 ***************************************************************************/
std::string GSPIModelDataSpace::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GSPIModelDataSpace ===");

        // Append attributes
        result.append("\n"+print_attributes());

        // Append parameter summary
        result.append("\n"+gammalib::parformat("Method"));
        result.append(m_method);
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
void GSPIModelDataSpace::init_members(void)
{
    // Initialise base class members
    m_instruments.clear();
    m_instruments.push_back("SPI");

    // Initialise members
    m_obs      = NULL;
    m_method.clear();
    m_index    = -1;
    m_map_size =  0;
    m_map      = NULL;
    m_parameters.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Model.
 ***************************************************************************/
void GSPIModelDataSpace::copy_members(const GSPIModelDataSpace& model)
{
    // Copy members
    m_obs        = model.m_obs;       // Copy pointer
    m_method     = model.m_method;
    m_index      = model.m_index;
    m_map_size   = model.m_map_size;
    m_parameters = model.m_parameters;

    // Copy parameter map
    if (m_map_size > 0) {
        m_map = new int[m_map_size];
        for (int i = 0; i < m_map_size; ++i) {
            m_map[i] = model.m_map[i];
        }
    }

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GSPIModelDataSpace::free_members(void)
{
    // Delete memory
    if (m_map != NULL) delete [] m_map;

    // Set pointers to free
    m_map = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set pointers
 *
 * Set pointers to all model parameters. The pointers are stored in a vector
 * that is member of the GModel base class.
 ***************************************************************************/
void GSPIModelDataSpace::set_pointers(void)
{
    // Clear parameter pointer(s)
    m_pars.clear();

    // Get number of parameters
    int npars = m_parameters.size();

    // Set pointers for all parameters
    for (int i = 0; i < npars; ++i) {
        m_pars.push_back(&(m_parameters[i]));
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Setup model
 *
 * @param[in] obs Observation.
 *
 * Setup the model for a given observation.
 ***************************************************************************/
void GSPIModelDataSpace::setup_model(const GObservation& obs) const
{
    // Continue only if observation pointer differs
    if (&obs != m_obs) {

        // Debug: signal that pointer differs
        #if defined(G_DEBUG_SETUP_MODEL)
        std::cout << "GSPIModelDataSpace::setup_model: ";
        std::cout << "observation pointer differs";
        std::cout << std::endl;
        #endif

        // Extract INTEGRAL/SPI observation
        const GSPIObservation* spi_obs = dynamic_cast<const GSPIObservation*>(&obs);
        if (spi_obs == NULL) {
            std::string cls = std::string(typeid(&obs).name());
            std::string msg = "Observation of type \""+cls+"\" is not an "
                              "INTEGRAL/SPI observation. Please specify an "
                              "INTEGRAL/SPI observation as argument.";
            throw GException::invalid_argument(G_SETUP_MODEL, msg);
        }

        // Store observation pointer
        m_obs = const_cast<GSPIObservation*>(spi_obs);

        // Extract INTEGRAL/SPI event cube. Continue only if the cube is
        // valid
        GSPIEventCube* cube = dynamic_cast<GSPIEventCube*>(m_obs->events());
        if (cube != NULL) {

            // Continue only if the model index is valid
            if (m_index >= 0 and m_index < cube->models()) {

                // Continue only if there is a method
                if (!m_method.empty()) {
                    const_cast<GSPIModelDataSpace*>(this)->setup_pars(cube);
                }

                // Debug: signal that index is not valid
                #if defined(G_DEBUG_SETUP_MODEL)
                else {
                    std::cout << "GSPIModelDataSpace::setup_model: ";
                    std::cout << "no method defined.";
                    std::cout << std::endl;
                }
                #endif

            } // endif: model index was valid

            // Debug: signal that index is not valid
            #if defined(G_DEBUG_SETUP_MODEL)
            else {
                std::cout << "GSPIModelDataSpace::setup_model: ";
                std::cout << "model index " << m_index << " is not in valid ";
                std::cout << "range.";
                std::cout << std::endl;
            }
            #endif

        } // endif: cube was valid

        // Debug: signal that event cube was invalid
        #if defined(G_DEBUG_SETUP_MODEL)
        else {
            std::cout << "GSPIModelDataSpace::setup_model: ";
            std::cout << "no valid event cube found in observation";
            std::cout << std::endl;
        }
        #endif

    } // endif: observation pointer differed

    // Return
    return;
}


/***********************************************************************//**
 * @brief Setup parameters
 *
 * @param[in] cube INTEGRAL/SPI event cube
 *
 * Setup the parameters for a given model.
 *
 * We assume that the model index is within the valid range of possible
 * indices and that a method was specified.
 ***************************************************************************/
void GSPIModelDataSpace::setup_pars(GSPIEventCube* cube)
{
    // Free parameter map
    free_members();

    // Get cube dimensions
    int n_pt  = cube->naxis(0);
    int n_det = cube->naxis(1);
    int n_eng = cube->naxis(2);

    // Initialise index vectors for all three dimensions
    std::vector<int> pt_indices(n_pt,   0);
    std::vector<int> det_indices(n_det, 0);
    std::vector<int> eng_indices(n_eng, 0);

    // Initialise name vectors for all three dimensions
    std::vector<std::string> pt_names;
    std::vector<std::string> det_names;
    std::vector<std::string> eng_names;

    // Determine pointing indices
    setup_pointing_indices(cube, &pt_indices, &pt_names);

    // Determine detector indices
    setup_detector_indices(cube, &det_indices, &det_names);

    // Determine energy indices
    setup_energy_indices(cube, &eng_indices, &eng_names);

    // Determine number of parameters in each dimension. The number of
    // parameters is given by the last index in the vector + 1. The check
    // for a non-zero index vector size if just for security, if the event
    // cube is set it should always be positive, yet for an empty event cube
    // it may be zero, and in this case the method will just do nothing.
    int npt   = (pt_indices.size()  > 0) ? pt_indices[pt_indices.size()-1]   + 1 : 0;
    int ndet  = (det_indices.size() > 0) ? det_indices[det_indices.size()-1] + 1 : 0;
    int neng  = (eng_indices.size() > 0) ? eng_indices[eng_indices.size()-1] + 1 : 0;
    int npars = npt * ndet * neng;

    // Debug: print number of parameters in each dimension
    #if defined(G_DEBUG_SETUP_MODEL)
    std::cout << "GSPIModelDataSpace::setup_pars: ";
    std::cout << "determined pointing, detector and energy indices:" << std::endl;
    std::cout << "- npt = " << npt << std::endl;
    std::cout << "- ndet = " << ndet << std::endl;
    std::cout << "- neng = " << neng << std::endl;
    std::cout << "- npars = " << npars << std::endl;
    #endif

    // Continue only if there are parameters in the model
    if (npars > 0) {

        // Allocate parameters vector
        std::vector<GModelPar> parameters;
        parameters.reserve(npars);

        // Allocate parameters. Energy parameter indices vary the most fastest
        // followed by detector parameter indices, follwed by pointing parameter
        // indices.
        for (int ipt = 0; ipt < npt; ++ipt) {
            for (int idet = 0; idet < ndet; ++idet) {
                for (int ieng = 0; ieng < neng; ++ieng) {

                    // Set energy name
                    std::string eng_name = (eng_names.size()  > 0) ? eng_names[ieng] : "";

                    // Build parameter name
                    std::string par_name;
                    if (m_name.empty()) {
                        par_name = "Scale";
                    }
                    else {
                        par_name = m_name;
                    }
                    if (eng_names.size() > 0) {
                        par_name += " " + eng_names[ieng];
                    }
                    if (det_names.size() > 0) {
                        par_name += " " + det_names[idet];
                    }
                    if (pt_names.size() > 0) {
                        par_name += " " + pt_names[ipt];
                    }

                    // If parameter name exists already in this class instance
                    // (because they were read using the read() method) then
                    // recover the parameter and push it into the vector of
                    // parameters
                    bool has_par = false;
                    for (int i = 0; i < m_parameters.size(); ++i) {
                        if (m_parameters[i].name() == par_name) {
                            parameters.push_back(m_parameters[i]);
                            has_par = true;
                            #if defined(G_DEBUG_SETUP_MODEL)
                            std::cout << "GSPIModelDataSpace::setup_pars: ";
                            std::cout << "copy old parameter \"";
                            std::cout << m_parameters[i].name() << "\"" << std::endl;
                            #endif
                            break;
                        }
                    }

                    // ... otherwise allocate new model parameter
                    if (!has_par) {
                        GModelPar par;
                        par.clear();
                        par.name(par_name);
                        par.free();
                        par.value(1.0);
                        par.scale(1.0);
                        par.gradient(0.0);
                        par.has_grad(true);
                        parameters.push_back(par);
                        #if defined(G_DEBUG_SETUP_MODEL)
                        std::cout << "GSPIModelDataSpace::setup_pars: ";
                        std::cout << "allocate new parameter \"";
                        std::cout << par_name << "\"" << std::endl;
                        #endif
                    }

                } // endfor: looped over energy parameters
            } // endfor: looped over detector parameters
        } // endfor: looped over pointing parameters

        // Allocate parameter map
        m_map_size = n_pt * n_det * n_eng;
        if (m_map_size > 0) {
            m_map = new int[m_map_size];
        }

        // Debug: print size of parameter map in each dimension
        #if defined(G_DEBUG_SETUP_MODEL)
        std::cout << "GSPIModelDataSpace::setup_pars: ";
        std::cout << "allocated parameter map:" << std::endl;
        std::cout << "- n_pt = " << n_pt << std::endl;
        std::cout << "- n_det = " << n_det << std::endl;
        std::cout << "- n_eng = " << n_eng << std::endl;
        std::cout << "- m_map_size = " << m_map_size << std::endl;
        #endif

        // Compute parameter map
        for (int i_pt = 0; i_pt < n_pt; ++i_pt) {
            for (int i_det = 0; i_det < n_det; ++i_det) {
                for (int i_eng = 0; i_eng < n_eng; ++i_eng) {

                    // Compute parameter index, respecting the scheme defined
                    // earlier: energy parameter indices vary the most fastest,
                    // followed by detector parameter indices, followed by
                    // pointing parameter indices
                    int index = pt_indices[i_pt]   * ndet * neng +
                                det_indices[i_det] * neng +
                                eng_indices[i_eng];

                    // Compute map index (same as in GSPIEventCube::read_dsp)
                    int i_map = (i_pt  * n_det + i_det) * n_eng + i_eng;

                    // Set map value
                    m_map[i_map] = index;

                    // Debug: print mapping
                    #if defined(G_DEBUG_SETUP_MODEL_MAP)
                    std::cout << "- [" << i_pt << "," << i_det << "," << i_eng << "]";
                    std::cout << " = " << index << std::endl;
                    #endif

                } // endfor: looped over energies
            } // endfor: looped over detectors
        } // endfor: looped over pointings

        // Replace model parameters
        m_parameters = parameters;

        // Set pointers
        set_pointers();

    } // endif: the model has parameters

    // Return
    return;
}


/***********************************************************************//**
 * @brief Setup pointing indices
 *
 * @param[in] cube Event cube.
 * @param[in,out] indices Vector of pointing indices.
 * @param[in,out] names Vector of pointing names.
 *
 * Setup vectors of pointing indices and pointing names. The length of the
 * vectors correspond to the number of pointings in the observation.
 ***************************************************************************/
void GSPIModelDataSpace::setup_pointing_indices(GSPIEventCube*            cube,
                                                std::vector<int>*         indices,
                                                std::vector<std::string>* names)
{
    // Convert method to lower case
    std::string method = gammalib::tolower(m_method);

    // Search for "point"
    if (gammalib::contains(method, "point")) {
        setup_point(cube, indices, names);
    }

    // Search for "orbit"
    else if (gammalib::contains(method, "orbit")) {
        setup_orbit(cube, indices, names);
    }

    // Search for "date"
    else if (gammalib::contains(method, "date")) {
        setup_date(cube, indices, names);
    }

    // Handle "gedfail"
    if (gammalib::contains(method, "gedfail")) {
        add_gedfail(cube, indices, names);
    }

    // Handle "gedanneal"
    if (gammalib::contains(method, "gedanneal")) {
        add_gedanneal(cube, indices, names);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Setup detector indices
 *
 * @param[in] cube Event cube.
 * @param[in,out] indices Vector of detector indices.
 * @param[in,out] names Vector of detector names.
 *
 * Setup vectors of detector indices and detector names. The length of the
 * vectors correspond to the number of detectors in the observation.
 ***************************************************************************/
void GSPIModelDataSpace::setup_detector_indices(GSPIEventCube*            cube,
                                                std::vector<int>*         indices,
                                                std::vector<std::string>* names)
{
    // Convert method to lower case
    std::string method = gammalib::tolower(m_method);

    // Search for "dete"
    if (gammalib::contains(method, "dete")) {
        setup_dete(cube, indices, names);
    }

    // Search for "evtclass"
    else if (gammalib::contains(method, "evtclass")) {
        setup_evtclass(cube, indices, names);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Setup energy indices
 *
 * @param[in] cube Event cube.
 * @param[in,out] indices Vector of energy indices.
 * @param[in,out] names Vector of energy names.
 *
 * Setup vectors of energy indices and energy names. The length of the
 * vectors correspond to the number of energies in the observation.
 ***************************************************************************/
void GSPIModelDataSpace::setup_energy_indices(GSPIEventCube*            cube,
                                              std::vector<int>*         indices,
                                              std::vector<std::string>* names)
{
    // Convert method to lower case
    std::string method = gammalib::tolower(m_method);

    // Search for "ebin"
    if (gammalib::contains(method, "ebin")) {
        setup_ebin(cube, indices, names);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Setup pointing indices and names for "point" method
 *
 * @param[in] cube Event cube.
 * @param[in,out] indices Vector of pointing indices.
 * @param[in,out] names Vector of pointing names.
 *
 * Setup vectors of pointing indices and names for "point" method.
 ***************************************************************************/
void GSPIModelDataSpace::setup_point(GSPIEventCube*            cube,
                                     std::vector<int>*         indices,
                                     std::vector<std::string>* names)
{
    // TODO: implement "point" method

    // Return
    return;
}


/***********************************************************************//**
 * @brief Setup pointing indices and names for "orbit" method
 *
 * @param[in] cube Event cube.
 * @param[in,out] indices Vector of pointing indices.
 * @param[in,out] names Vector of pointing names.
 *
 * Setup vectors of pointing indices and names for "orbit" method.
 ***************************************************************************/
void GSPIModelDataSpace::setup_orbit(GSPIEventCube*            cube,
                                     std::vector<int>*         indices,
                                     std::vector<std::string>* names)
{
    // Allocate orbit strings
    std::vector<std::string> orbits;

    // Get number of pointings
    int npt = indices->size();

    // Loop over all pointings
    for (int ipt = 0; ipt < npt; ++ipt) {

        // Initialise index
        int index = -1;

        // Get orbit
        std::string orbit = cube->ptid(ipt).substr(0,4);

        // If orbit is in orbits then store the location of the string
        // in the vector as parameter index.
        for (int i = 0; i < orbits.size(); ++i) {
            if (orbit == orbits[i]) {
                index = i;
                break;
            }
        }

        // If orbit is not yet in vector then
        if (index == -1) {
            index = orbits.size();
            orbits.push_back(orbit);
            names->push_back("Rev" + orbit);
        }

        // Store index
        (*indices)[ipt] = index;

    } // endfor: looped over pointings

    // Return
    return;
}


/***********************************************************************//**
 * @brief Setup pointing indices and names for "date" method
 *
 * @param[in] cube Event cube.
 * @param[in,out] indices Vector of pointing indices.
 * @param[in,out] names Vector of pointing names.
 *
 * Setup vectors of pointing indices and names for "date" method.
 ***************************************************************************/
void GSPIModelDataSpace::setup_date(GSPIEventCube*            cube,
                                    std::vector<int>*         indices,
                                    std::vector<std::string>* names)
{
    // TODO: implement "date" method

    // Return
    return;
}


/***********************************************************************//**
 * @brief Modify pointing indices and names for "gedfail" method
 *
 * @param[in] cube Event cube.
 * @param[in,out] indices Vector of pointing indices.
 * @param[in,out] names Vector of pointing names.
 *
 * Inserts additional model parameters for each pointing where a detector
 * failure occured.
 ***************************************************************************/
void GSPIModelDataSpace::add_gedfail(GSPIEventCube*            cube,
                                     std::vector<int>*         indices,
                                     std::vector<std::string>* names)
{
    // TODO: implement "gedfail" method

    // Return
    return;
}


/***********************************************************************//**
 * @brief Modify pointing indices and names for "gedanneal" method
 *
 * @param[in] cube Event cube.
 * @param[in,out] indices Vector of pointing indices.
 * @param[in,out] names Vector of pointing names.
 *
 * Inserts additional model parameters for each pointing where a detector
 * annealing occured.
 ***************************************************************************/
void GSPIModelDataSpace::add_gedanneal(GSPIEventCube*            cube,
                                       std::vector<int>*         indices,
                                       std::vector<std::string>* names)
{
    // TODO: implement "gedanneal" method

    // Return
    return;
}


/***********************************************************************//**
 * @brief Setup detector indices and names for "dete" method
 *
 * @param[in] cube Event cube.
 * @param[in,out] indices Vector of detector indices.
 * @param[in,out] names Vector of detector names.
 *
 * Setup vectors of detector indices and names for "dete" method.
 ***************************************************************************/
void GSPIModelDataSpace::setup_dete(GSPIEventCube*            cube,
                                    std::vector<int>*         indices,
                                    std::vector<std::string>* names)
{
    // Get number of detectors
    int ndet = indices->size();

    // Setup index vector
    for (int idet = 0; idet < ndet; ++idet) {
        (*indices)[idet] = idet;
        names->push_back("D" + gammalib::str(idet, "%3.3d"));
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Setup detector indices and names for "evtclass" method
 *
 * @param[in] cube Event cube.
 * @param[in,out] indices Vector of detector indices.
 * @param[in,out] names Vector of detector names.
 *
 * Setup vectors of detector indices and names for "evtclass" method.
 ***************************************************************************/
void GSPIModelDataSpace::setup_evtclass(GSPIEventCube*            cube,
                                        std::vector<int>*         indices,
                                        std::vector<std::string>* names)
{
    // TODO: implement "evtclass" method

    // Return
    return;
}


/***********************************************************************//**
 * @brief Setup energy indices and names for "ebin" method
 *
 * @param[in] cube Event cube.
 * @param[in,out] indices Vector of energy indices.
 * @param[in,out] names Vector of energy names.
 *
 * Setup vectors of energy indices and names for "ebin" method.
 ***************************************************************************/
void GSPIModelDataSpace::setup_ebin(GSPIEventCube*            cube,
                                    std::vector<int>*         indices,
                                    std::vector<std::string>* names)
{
    // Get number of detectors
    int neng = indices->size();

    // Setup index vector
    for (int ieng = 0; ieng < neng; ++ieng) {
        (*indices)[ieng] = ieng;
        names->push_back("E" + gammalib::str(ieng, "%3.3d"));
    }

    // Return
    return;
}
