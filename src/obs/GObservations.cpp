/***************************************************************************
 *                GObservations.cpp - Observation container class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2013 by Juergen Knoedlseder                         *
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
 * @file GObservations.cpp
 * @brief Observations container class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GObservations.hpp"
#include "GObservationRegistry.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_AT                                        "GObservations::at(int&)"
#define G_SET                       "GObservations::set(int&, GObservation&)"
#define G_APPEND                       "GObservations::append(GObservation&)"
#define G_INSERT                 "GObservations::insert(int&, GObservation&)"
#define G_REMOVE                                "GObservations::remove(int&)"
#define G_EXTEND                      "GObservations::extend(GObservations&)"
#define G_READ                                   "GObservations::read(GXml&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GObservations::GObservations(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param obs Observation container.
 ***************************************************************************/
GObservations::GObservations(const GObservations& obs)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(obs);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load constructor
 *
 * @param[in] filename XML filename.
 *
 * Construct the observation container by loading all relevant information
 * from a XML file. Please refer to the read() method for more information
 * about the structure of the XML file.
 ***************************************************************************/
GObservations::GObservations(const std::string& filename)
{
    // Initialise members
    init_members();

    // Load XML file
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GObservations::~GObservations(void)
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
 * @param[in] obs Observation container.
 * @return Observation container.
 ***************************************************************************/
GObservations& GObservations::operator=(const GObservations& obs)
{
    // Execute only if object is not identical
    if (this != &obs) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(obs);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear observations
 ***************************************************************************/
void GObservations::clear(void)
{
    // Free members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone observations
 *
 * @return Pointer to deep copy of observation container.
 ***************************************************************************/
GObservations* GObservations::clone(void) const
{
    // Clone observations
    return new GObservations(*this);
}


/***********************************************************************//**
 * @brief Return pointer to observation
 *
 * @param[in] index Observation index [0,...,size()-1].
 * @return Observation.
 *
 * @exception GException::out_of_range
 *            Operation index is out of range.
 *
 * Returns a pointer to the observation with specified @p index.
 ***************************************************************************/
GObservation* GObservations::at(const int& index)
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, index, 0, size()-1);
    }

    // Return pointer
    return m_obs[index];
}


/***********************************************************************//**
 * @brief Return pointer to observation (const version)
 *
 * @param[in] index Observation index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Operation index is out of range.
 *
 * Returns a const pointer to the observation with specified @p index.
 ***************************************************************************/
const GObservation* GObservations::at(const int& index) const
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, index, 0, size()-1);
    }

    // Return pointer
    return m_obs[index];
}


/***********************************************************************//**
 * @brief Set observation in container
 *
 * @param[in] index Observation index [0,...,size()-1].
 * @param[in] obs Observation.
 * @return Pointer to deep copy of observation.
 *
 * @exception GException::out_of_range
 *            Observation index is out of range.
 * @exception GException::invalid_value
 *            Observation with same instrument and identifier already
 *            exists in container.
 *
 * Set a deep copy and observation @p obs at the specified @p index in the
 * container.
 ***************************************************************************/
GObservation* GObservations::set(const int& index, const GObservation& obs)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_SET, index, 0, size()-1);
    }
    #endif

    // Raise an exception if an observation with specified instrument and
    // identifier already exists
    int inx = get_index(obs.instrument(), obs.id());
    if (inx != -1 && inx != index) {
        std::string msg =
            "Attempt to set \""+obs.instrument()+"\" observation with"
            " identifier \""+obs.id()+"\" in observation container at"
            " index "+gammalib::str(index)+", but an observation with the"
            " same attributes exists already at index "+gammalib::str(inx)+
            " in the container.\n"
            "Every observation for a given instrument in the observation"
            " container needs a unique identifier.";
        throw GException::invalid_value(G_SET, msg);
    }

    // Delete existing observation
    if (m_obs[index] != NULL) delete m_obs[index];

    // Store pointer to a deep copy of the observation
    m_obs[index] = obs.clone();

    // Return pointer
    return m_obs[index];
}


/***********************************************************************//**
 * @brief Append observation to container
 *
 * @param[in] obs Observation.
 * @return Pointer to deep copy of observation.
 *
 * @exception GException::invalid_value
 *            Observation with same instrument and identifier already
 *            exists in container.
 *
 * Appends a deep copy of an observation to the container.
 ***************************************************************************/
GObservation* GObservations::append(const GObservation& obs)
{
    // Raise an exception if an observation with specified instrument and
    // identifier already exists
    int inx = get_index(obs.instrument(), obs.id());
    if (inx != -1) {
        std::string msg =
            "Attempt to append \""+obs.instrument()+"\" observation with"
            " identifier \""+obs.id()+"\" to observation container, but an"
            " observation with the same attributes exists already at"
            " index "+gammalib::str(inx)+" in the container.\n"
            "Every observation for a given instrument in the observation"
            " container needs a unique identifier.";
        throw GException::invalid_value(G_APPEND, msg);
    }

    // Clone observation
    GObservation* ptr = obs.clone();
    
    // Append to list
    m_obs.push_back(ptr);

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Insert observation into container
 *
 * @param[in] index Observation index [0,...,size()-1].
 * @param[in] obs Observation.
 * @return Pointer to deep copy of observation.
 *
 * @exception GException::out_of_range
 *            Observation index is out of range.
 * @exception GException::invalid_value
 *            Observation with same instrument and identifier already
 *            exists in container.
 *
 * Inserts a deep copy of an observation into the container before the
 * observation with the specified @p index.
 ***************************************************************************/
GObservation* GObservations::insert(const int& index, const GObservation& obs)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (isempty()) {
        if (index > 0) {
            throw GException::out_of_range(G_INSERT, index, 0, size()-1);
        }
    }
    else {
        if (index < 0 || index >= size()) {
            throw GException::out_of_range(G_INSERT, index, 0, size()-1);
        }
    }
    #endif

    // Raise an exception if an observation with specified instrument and
    // identifier already exists
    int inx = get_index(obs.instrument(), obs.id());
    if (inx != -1 && inx != index) {
        std::string msg =
            "Attempt to insert \""+obs.instrument()+"\" observation with"
            " identifier \""+obs.id()+"\" in observation container before"
            " index "+gammalib::str(index)+", but an observation with the"
            " same attributes exists already at index "+gammalib::str(inx)+
            " in the container.\n"
            "Every observation for a given instrument in the observation"
            " container needs a unique identifier.";
        throw GException::invalid_value(G_INSERT, msg);
    }

    // Clone observation
    GObservation* ptr = obs.clone();

    // Clone observation and insert into list
    m_obs.insert(m_obs.begin()+index, ptr);

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Remove observation from container
 *
 * @param[in] index Observation index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Observation index is out of range.
 *
 * Removes observation of specified @p index from the container.
 ***************************************************************************/
void GObservations::remove(const int& index)
{
    // Compile option: If index is outside boundary then raise exception
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_REMOVE, index, 0, size()-1);
    }
    #endif

    // Erase observation component from container
    m_obs.erase(m_obs.begin() + index);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append observations from observation container
 *
 * @param[in] obs Observations.
 *
 * @exception GException::invalid_value
 *            Observation with same instrument and identifier already
 *            exists in container.
 *
 * Appends deep copies of observations to the container.
 ***************************************************************************/
void GObservations::extend(const GObservations& obs)
{
    // Do nothing if observation container is empty
    if (!obs.isempty()) {

        // Get size. Note that we extract the size first to avoid an
        // endless loop that arises when a container is appended to
        // itself.
        int num = obs.size();

        // Reserve enough space
        reserve(size() + num);

        // Loop over all observations, clone them and append them to the
        // list
        for (int i = 0; i < num; ++i) {

            // Raise an exception if an observation with specified
            // instrument and identifier already exists
            int inx = get_index(obs[i]->instrument(), obs[i]->id());
            if (inx != -1) {
                std::string msg =
                    "Attempt to append \""+obs[i]->instrument()+"\""
                    " observation with identifier \""+obs[i]->id()+"\""
                    " to observation container, but an observation with the"
                    " same attributes exists already at index "+
                    gammalib::str(inx)+" in the container.\n"
                    "Every observation for a given instrument in the"
                    " observation container needs a unique identifier.";
                throw GException::invalid_value(G_EXTEND, msg);
            }

            // Append observation to container
            m_obs.push_back(obs[i]->clone());
        }

    } // endif: observation container was not empty
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Signals if observation exists
 *
 * @param[in] instrument Instrument.
 * @param[in] id Observation identifier.
 * @return True if observation with specified @p instrument and identifier
 *         @p id exists.
 *
 * Searches all observations for a match with the specified @p instrument
 * and @p identifier. If the specified attributes have been found, true is
 * returned.
 ***************************************************************************/
bool GObservations::contains(const std::string& instrument,
                             const std::string& id) const
{
    // Get observation index
    int index = get_index(instrument, id);

    // Return
    return (index != -1);
}


/***********************************************************************//**
 * @brief Load observations from XML file
 *
 * @param[in] filename XML filename.
 *
 * Loads observation from a XML file into the container. Please refer to the
 * read() method for more information about the structure of the XML file.
 ***************************************************************************/
void GObservations::load(const std::string& filename)
{
    // Clear any existing observations
    clear();

    // Load XML document
    GXml xml(filename);

    // Read observations from XML document
    read(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save observations into XML file
 *
 * @param[in] filename XML filename.
 *
 * Saves observations into a XML file. Please refer to the read() method for
 * more information about the structure of the XML file.
 ***************************************************************************/
void GObservations::save(const std::string& filename) const
{
    // Declare empty XML document
    GXml xml;

    // Write observations into XML file
    write(xml);

    // Save XML document
    xml.save(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read observations from XML document
 *
 * @param[in] xml XML document.
 *
 * @exception GException::invalid_instrument
 *            Invalid instrument encountered in XML file.
 *
 * Reads observations from the first observation list that is found in the
 * XML document. The decoding of the instrument specific observation
 * definition is done within the observation's GObservation::read() method.
 * The following file structure is expected:
 *
 *     <observation_list title="observation library">
 *       <observation name="..." id="..." instrument="...">
 *         ...
 *       </observation>
 *       <observation name="..." id="..." instrument="...">
 *         ...
 *       </observation>
 *       ...
 *     </observation_list>
 *
 * The @p name and @p id attributes allow for a unique identification of an
 * observation within the observation container. The @p instrument
 * attributes specifies the instrument to which the observation applies.
 * This attribute will be used to allocate the appropriate instrument
 * specific GObservation class variant by using the GObservationRegistry
 * class.
 *
 * The structure within the @p observation tag is defined by the instrument
 * specific GObservation class.
 *
 * @todo Observation names and IDs are not verified so far for uniqueness.
 *       This would be required to achieve an unambiguous update of parameters
 *       in an already existing XML file when using the write method.
 ***************************************************************************/
void GObservations::read(const GXml& xml)
{
    // Get pointer on observation library
    const GXmlElement* lib = xml.element("observation_list", 0);

    // Loop over all observations
    int n = lib->elements("observation");
    for (int i = 0; i < n; ++i) {

        // Get pointer on observation
        const GXmlElement* obs = lib->element("observation", i);

        // Get attributes
        std::string name       = obs->attribute("name");
        std::string id         = obs->attribute("id");
        std::string instrument = obs->attribute("instrument");

        // Allocate instrument specific observation
        GObservationRegistry registry;
        GObservation*        ptr = registry.alloc(instrument);

        // If observation is valid then read its definition from XML file
        if (ptr != NULL) {

            // Read definition
            ptr->read(*obs);

            // Set attributes
            ptr->name(name);
            ptr->id(id);

        } // endif: observation was valid

        // ... otherwise throw an exception
        else {
            throw GException::invalid_instrument(G_READ, instrument);
        }

        // Append observation to container
        append(*ptr);

        // Free observation (the append method made a deep copy)
        delete ptr;

    } // endfor: looped over all observations

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write observations into XML document
 *
 * @param[in] xml XML document.
 *
 * Write observations into the first observation library that is found in the
 * XML document. In case that no observation library exists, one is added to
 * the document. Please refer to the read() method for more information
 * about the structure of the XML document.
 ***************************************************************************/
void GObservations::write(GXml& xml) const
{
    // If there is no observation library then append one
    if (xml.elements("observation_list") == 0) {
        xml.append(GXmlElement("observation_list title=\"observation list\""));
    }

    // Get pointer on observation library
    GXmlElement* lib = xml.element("observation_list", 0);

    // Write all observations into library
    for (int i = 0; i < size(); ++i) {
    
        // Initialise pointer on observation
        GXmlElement* obs = NULL;

        // Search corresponding observation
        int n = xml.elements("observation");
        for (int k = 0; k < n; ++k) {
            GXmlElement* element = xml.element("observation", k);
            if (element->attribute("name")       == m_obs[i]->name() &&
                element->attribute("id")         == m_obs[i]->id()   &&
                element->attribute("instrument") == m_obs[i]->instrument()) {
                obs = element;
                break;
            }
        }

        // If no observation with corresponding name, ID and instrument was
        // found then append one now
        if (obs == NULL) {
            obs = lib->append("observation");
            obs->attribute("name", m_obs[i]->name());
            obs->attribute("id", m_obs[i]->id());
            obs->attribute("instrument", m_obs[i]->instrument());
        }
    
        // Write now observation
        m_obs[i]->write(*obs);

    } // endfor: looped over all observaitons

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load models from XML file
 *
 * @param[in] filename XML filename.
 *
 * Loads the models from a XML file. Please refer to the GModels::read()
 * method for more information about the expected structure of the XML
 * file. 
 ***************************************************************************/
void GObservations::models(const std::string& filename)
{
    // Load models
    m_models.load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Optimize model parameters using optimizer
 *
 * @param[in] opt Optimizer.
 *
 * Optimizes the free parameters of the models by using the optimizer
 * that has been provided by the @p opt argument.
 ***************************************************************************/
void GObservations::optimize(GOptimizer& opt)
{
    // Extract optimizer parameter container from model container
    GOptimizerPars pars = m_models.pars();

    // Optimize model parameters
    opt.optimize(m_fct, pars);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print observation list information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing observation list information
 ***************************************************************************/
std::string GObservations::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GObservations ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of observations"));
        result.append(gammalib::str(size()));
        result.append("\n"+gammalib::parformat("Number of predicted events"));
        result.append(gammalib::str(npred()));

        // NORMAL: Append observations
        if (chatter >= NORMAL) {
            for (int i = 0; i < size(); ++i) {
                result.append("\n");
                result.append((*this)[i]->print());
            }
        } // endif: chatter was normal

        // EXPLICIT: Append models
        if (chatter >= EXPLICIT) {
            result.append("\n"+m_models.print());
        } // endif: chatter was explicit

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                              Private methods                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GObservations::init_members(void)
{
    // Initialise members
    m_obs.clear();
    m_models.clear();
    m_fct.set(this);  //!< Makes sure that optimizer points to this instance

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] obs Observation container.
 *
 * Copies all class members. Deep copies are created for the observations,
 * which is what is needed to have a real physical copy of the members.
 ***************************************************************************/
void GObservations::copy_members(const GObservations& obs)
{
    // Copy attributes. WARNING: The member m_fct SHALL not be copied to not
    // corrupt its m_this pointer which should always point to the proper
    // observation. See note in init_members().
    m_models = obs.m_models;

    // Copy observations
    m_obs.clear();
    for (int i = 0; i < obs.m_obs.size(); ++i) {
        m_obs.push_back((obs.m_obs[i]->clone()));
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 *
 * Deallocates all observations. Since container classes that hold pointers
 * need to handle the proper deallocation of memory, we loop here over all
 * pointers and make sure that we deallocate all observations in the
 * container.
 ***************************************************************************/
void GObservations::free_members(void)
{
    // Free observations
    for (int i = 0; i < m_obs.size(); ++i) {
        delete m_obs[i];
        m_obs[i] = NULL;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return observation index by instrument and identifier
 *
 * @param[in] instrument Instrument.
 * @param[in] id Observation identifier.
 * @return Observation index (-1 of instrument and observation identifier
 *         has not been found)
 *
 * Returns observation index based on the specified @p instrument and
 * observation identifier @p id. If no observation with the specified
 * attributes has been found, the method returns -1.
 ***************************************************************************/
int GObservations::get_index(const std::string& instrument,
                             const std::string& id) const
{
    // Initialise index
    int index = -1;

    // Search observation with specified instrument and id
    for (int i = 0; i < size(); ++i) {
        if ((m_obs[i]->instrument() == instrument) &&
            (m_obs[i]->id()         == id)) {
            index = i;
            break;
        }
    }

    // Return index
    return index;
}
