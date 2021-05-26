/***************************************************************************
 *              GModel.cpp - Abstract virtual model base class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2021 by Juergen Knoedlseder                         *
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
 * @file GModel.cpp
 * @brief Abstract model base class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GModel.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ACCESS                           "GModel::operator[](std::string&)"
#define G_AT                                               "GModel::at(int&)"
#define G_SCALE                              "GModelPar& GModel::scale(int&)"
#define G_WRITE_SCALES                   "GModel::write_scales(GXmlElement&)"

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
GModel::GModel(void)
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
 * Construct model from XML element. The method extracts all model attributes
 * from the XML file (see the read_attributes() method for more information
 * about the supported attributes).
 ***************************************************************************/
GModel::GModel(const GXmlElement& xml)
{
    // Initialise members
    init_members();

    // Read attributes
    read_attributes(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Model.
 ***************************************************************************/
GModel::GModel(const GModel& model)
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
GModel::~GModel(void)
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
 * @param[in] model Model.
 * @return Model.
 ***************************************************************************/
GModel& GModel::operator=(const GModel& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(model);

    } // endif: object was not identical

    // Return
    return *this;
}


/***********************************************************************//**
 * @brief Returns reference to model parameter by name
 *
 * @param[in] name Parameter name.
 * @return Reference to model parameter.
 *
 * @exception GException::invalid_argument
 *            Parameter with specified name not found.
 *
 * Returns a reference to the model parameter of the specified @p name.
 * Throws an exception if no parameter with @p name is found.
 ***************************************************************************/
GModelPar& GModel::operator[](const std::string& name)
{
    // Get parameter index
    int index = 0;
    for (; index < size(); ++index) {
        if (m_pars[index]->name() == name) {
            break;
        }
    }

    // Throw exception if parameter name was not found
    if (index >= size()) {
        std::string msg = "Model parameter \""+name+"\" not found in model. "
                          "Please specify a valid model parameter name.";
        throw GException::invalid_argument(G_ACCESS, msg);
    }

    // Return reference
    return *(m_pars[index]);
}


/***********************************************************************//**
 * @brief Returns reference to model parameter by name const version)
 *
 * @param[in] name Parameter name.
 * @return Reference to model parameter.
 *
 * @exception GException::invalid_argument
 *            Parameter with specified name not found.
 *
 * Returns a const reference to the model parameter of the specified
 * @p name. Throws an exception if no parameter with @p name is found.
 ***************************************************************************/
const GModelPar& GModel::operator[](const std::string& name) const
{
    // Get parameter index
    int index = 0;
    for (; index < size(); ++index) {
        if (m_pars[index]->name() == name) {
            break;
        }
    }

    // Throw exception if parameter name was not found
    if (index >= size()) {
        std::string msg = "Model parameter \""+name+"\" not found in model. "
                          "Please specify a valid model parameter name.";
        throw GException::invalid_argument(G_ACCESS, msg);
    }

    // Return reference
    return *(m_pars[index]);
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Returns reference to model parameter by index
 *
 * @param[in] index Parameter index [0,...,size()[.
 * @return Reference to model parameter.
 *
 * @exception GException::out_of_range
 *            Parameter index is out of range.
 *
 * Returns a reference to the model parameter of the specified @p index.
 * Throws an exception if @p index is not valid.
 ***************************************************************************/
GModelPar& GModel::at(const int& index)
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Parameter index", index, size());
    }

    // Return reference
    return *(m_pars[index]);
}


/***********************************************************************//**
 * @brief Returns reference to model parameter by index (const version)
 *
 * @param[in] index Parameter index [0,...,size()[.
 * @return Const reference to model parameter.
 *
 * @exception GException::out_of_range
 *            Parameter index is out of range.
 *
 * Returns a const reference to the model parameter of the specified
 * @p index. Throws an exception if @p index is not valid.
 ***************************************************************************/
const GModelPar& GModel::at(const int& index) const
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Parameter index", index, size());
    }

    // Return reference
    return *(m_pars[index]);
}


/***********************************************************************//**
 * @brief Checks if parameter name exists
 *
 * @param[in] name Parameter name.
 * @return True if parameter with specified @p name exists.
 *
 * Searches all parameter names for a match with the specified @p name. If
 * the specified name has been found, true is returned.
 ***************************************************************************/
bool GModel::has_par(const std::string& name) const
{
    // Default found flag to false
    bool found = false;

    // Search for parameter name
    for (int i = 0; i < size(); ++i) {
        if (m_pars[i]->name() == name) {
            found = true;
            break;
        }
    }

    // Return
    return found;
}


/***********************************************************************//**
 * @brief Returns instruments to which model applies
 *
 * @return Instruments.
 *
 * Returns a comma separated list of instruments to which model applies. If
 * no instrument exists then an empty string is returned.
 ***************************************************************************/
std::string GModel::instruments(void) const
{
    // Initialise string
    std::string result;

    // Attach all instruments
    for (int i = 0; i < m_instruments.size(); ++i) {
        if (i > 0) {
            result += ",";
        }
        result += m_instruments[i];
    }

    // Return
    return result;
}


/***********************************************************************//**
 * @brief Set instruments to which model applies
 *
 * @param[in] instruments String of instruments.
 *
 * Sets the instruments to which the model applies from a comma separated
 * list of strings. The instrument names are case sensitive.
 *
 * If the @p instrument string is empty, the model is considered to apply to
 * all instruments.
 ***************************************************************************/
void GModel::instruments(const std::string& instruments)
{
    // Clear instruments vector
    m_instruments.clear();

    // Extract instruments
    std::vector<std::string> inst = gammalib::split(instruments, ",");

    // Attach all instruments
    for (int i = 0; i < inst.size(); ++i) {
        m_instruments.push_back(gammalib::strip_whitespace(inst[i]));
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns reference to scale parameter by index
 *
 * @param[in] index Scale parameter index [0,...,scales()[.
 * @return Reference to scale parameter.
 *
 * @exception GException::out_of_range
 *            Scale parameter index is out of range.
 *
 * Returns a reference to the scale parameter of the specified @p index.
 * Throws an exception if @p index is not valid.
 ***************************************************************************/
GModelPar& GModel::scale(const int& index)
{
    // Raise exception if index is out of range
    if (index < 0 || index >= scales()) {
        throw GException::out_of_range(G_SCALE, "Scale parameter index",
                                       index, scales());
    }

    // Return reference
    return (m_scales[index]);
}


/***********************************************************************//**
 * @brief Returns reference to scale parameter by index (const version)
 *
 * @param[in] index Scale parameter index [0,...,scales()[.
 * @return Reference to scale parameter.
 *
 * @exception GException::out_of_range
 *            Scale parameter index is out of range.
 *
 * Returns a reference to the scale parameter of the specified @p index.
 * Throws an exception if @p index is not valid.
 ***************************************************************************/
const GModelPar& GModel::scale(const int& index) const
{
    // Raise exception if index is out of range
    if (index < 0 || index >= scales()) {
        throw GException::out_of_range(G_SCALE, "Scale parameter index",
                                       index, scales());
    }

    // Return reference
    return (m_scales[index]);
}


/***********************************************************************//**
 * @brief Returns model scale factor for a given instrument
 *
 * @param[in] instrument Instrument.
 *
 * Returns the model scale factor for a given @p instrument. The search is
 * case sensitive.
 *
 * If the @p instrument is not found, the method returns a scale factor of
 * unity.
 ***************************************************************************/
GModelPar GModel::scale(const std::string& instrument) const
{
    // Initialise unit scale factor
    GModelPar scale;
    scale.value(1.0);
    scale.name(instrument);
    scale.fix();

    // Search for instrument and recover scale factor if the instrument
    // has been found.
    for (int i = 0; i < m_scales.size(); ++i) {
        if (m_scales[i].name() == instrument) {
            scale = m_scales[i];
            break;
        }
    }

    // Return scale factor
    return scale;
}


/***********************************************************************//**
 * @brief Set model scale factor for a given instrument
 *
 * @param[in] par Model parameter for scaling.
 *
 * Sets the model parameter for a given instrument. The instrument name is
 * case sensitive, but any leading to trailing white space will be stripped.
 *
 * If the instrument is not yet defined it will be appended to the list of
 * instruments.
 ***************************************************************************/
void GModel::scale(const GModelPar& par)
{
    // String leading and trailing while space
    std::string instrument = gammalib::strip_whitespace(par.name());

    // Search for instrument and copy the model parameter if the instrument
    // has been found. Make sure that the instrument name is in upper case.
    bool found = false;
    for (int i = 0; i < m_scales.size(); ++i) {
        if (m_scales[i].name() == instrument) {
            found       = true;
            m_scales[i] = par;
            m_scales[i].name(instrument);
            m_scales[i].has_grad(true);
            break;
        }
    }

    // If instrument has not been found then append it now to the list
    // of instruments.
    if (!found) {

        // Push scale parameter in list
        m_scales.push_back(par);

        // Get index of last scale parameter
        int i = m_scales.size()-1;

        // Set instrument name and signal availability of gradient
        m_scales[i].name(instrument);
        m_scales[i].has_grad(true);

        // Push scale parameter on parameter stack
        m_pars.push_back(&m_scales[i]);

    } // endif: new scale parameter

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns observation identifiers to which model applies
 *
 * @return Observation identifiers.
 *
 * Returns a comma separated list of observation identifiers to which model
 * applies. If no observation identifier exists then an empty string is
 * returned.
 ***************************************************************************/
std::string GModel::ids(void) const
{
    // Initialise string
    std::string result;

    // Attach all observation identifiers
    for (int i = 0; i < m_ids.size(); ++i) {
        if (i > 0) {
            result += ",";
        }
        result += m_ids[i];
    }

    // Return
    return result;
}


/***********************************************************************//**
 * @brief Set observation identifiers to which model applies
 *
 * @param[in] ids String of observation identifiers.
 *
 * Sets the observation identifiers to which the model applies from a comma
 * separated list of strings. The observation identifiers are case sensitive,
 * but any leading and trailing whitespace will be stripped.
 *
 * If the observation identifier string is empty, the model is considered to
 * apply to all observation identifiers.
 ***************************************************************************/
void GModel::ids(const std::string& ids)
{
    // Clear observation identifier vector
    m_ids.clear();

    // Extract observation identifiers
    std::vector<std::string> id = gammalib::split(ids, ",");

    // Attach all observation identifiers
    for (int i = 0; i < id.size(); ++i) {
        m_ids.push_back(gammalib::strip_whitespace(id[i]));
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Verifies if model is valid for a given instrument and identifier
 *
 * @param[in] instrument Instrument name.
 * @param[in] id Observation identifier.
 * @return Validity flag.
 *
 * Checks if specified instrument name and observation identifier is in list
 * of applicable instruments and identifiers. The check is case sensitive.
 *
 * If the list of applicable instruments is empty, the model applies to all
 * possible instruments. If the list of applicable observation identifiers
 * is empty, the model applies to all identifiers.
 *
 * If an empty string is provided as @p instrument parameter, the check for
 * instrument validity will be skipped. Similarily, if an empty string is
 * provided as @p id parameter, the identifier check will be skipped. This
 * allows for example to check for models of a given identifier whatever the
 * instrument, or for models for a given instrument, whatever the identifier.
 ***************************************************************************/
bool GModel::is_valid(const std::string& instrument,
                      const std::string& id) const
{
    // Initialise validity
    bool valid = true;

    // Check if model applies to instrument
    if (!m_instruments.empty() && !instrument.empty()) {

        // Initialise validity flag
        valid = false;

        // Check if instrument is in list
        for (int i = 0; i < m_instruments.size(); ++i) {
            if (instrument == m_instruments[i]) {
                valid = true;
                break;
            }
        }

    }

    // Check if model applies to observation identifier
    if (valid && !m_ids.empty() && !id.empty()) {

        // Initialise validity flag
        valid = false;

        // Check if name is in list
        for (int i = 0; i < m_ids.size(); ++i) {
            if (id == m_ids[i]) {
                valid = true;
                break;
            }
        }

    }

    // Return
    return valid;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModel::init_members(void)
{
    // Initialise members
    m_name.clear();
    m_instruments.clear();
    m_scales.clear();
    m_ids.clear();
    m_pars.clear();
    m_associations.clear();
    m_ts           = 0.0;
    m_has_ts       = false;
    m_has_tscalc   = false;
    m_tscalc       = false;
    m_has_eval_inx = false;
    m_eval_inx.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Model.
 ***************************************************************************/
void GModel::copy_members(const GModel& model)
{
    // Copy members
    m_name         = model.m_name;
    m_instruments  = model.m_instruments;
    m_scales       = model.m_scales;
    m_ids          = model.m_ids;
    m_pars         = model.m_pars;
    m_associations = model.m_associations;
    m_ts           = model.m_ts;
    m_has_ts       = model.m_has_ts;
    m_has_tscalc   = model.m_has_tscalc;
    m_tscalc       = model.m_tscalc;
    m_has_eval_inx = model.m_has_eval_inx;
    m_eval_inx     = model.m_eval_inx;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModel::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Read model attributes
 *
 * @param[in] xml XML element.
 ***************************************************************************/
void GModel::read_attributes(const GXmlElement& xml)
{
    // Set model name
    name(xml.attribute("name"));

    // Set instruments
    instruments(xml.attribute("instrument"));

    // Set observation identifiers
    ids(xml.attribute("id"));

    // Set model TS
    if (xml.has_attribute("ts")) {
        std::string ts = xml.attribute("ts");
        this->ts(gammalib::todouble(ts));
    }

    // Set TS computation flag
    if (xml.has_attribute("tscalc")) {
        bool tscalc = (xml.attribute("tscalc") == "1") ? true : false;
        this->tscalc(tscalc);
    }

    // Read instrument scales
    read_scales(xml);

    // Read model associations
    m_associations.read(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model attributes
 *
 * @param[in] xml XML element.
 ***************************************************************************/
void GModel::write_attributes(GXmlElement& xml) const
{
    // Set model name
    xml.attribute("name", name());

    // Set model type
    xml.attribute("type", type());

    // Set instruments
    std::string instruments = this->instruments();
    if (instruments.length() > 0) {
        xml.attribute("instrument", instruments);
    }

    // Set observation identifiers
    std::string identifiers = ids();
    if (identifiers.length() > 0) {
        xml.attribute("id", identifiers);
    }

    // If available, set "ts" attribute
    if (m_has_ts) {
        xml.attribute("ts", gammalib::str(ts(), 3));
    }

    // If available, set "tscalc" attribute
    if (m_has_tscalc) {
        std::string ts_calc = tscalc() ? "1" : "0";
        xml.attribute("tscalc", ts_calc);
    }

    // Write instrument scales
    write_scales(xml);

    // Write model associations
    m_associations.write(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print model attributes
 *
 * @return Returns string with model attributes.
 ***************************************************************************/
std::string GModel::print_attributes(void) const
{
    // Initialise result string
    std::string result;

    // Append model name
    result.append(gammalib::parformat("Name")+name());

    // Append instruments
    result.append("\n"+gammalib::parformat("Instruments"));
    if (!m_instruments.empty()) {
        for (int i = 0; i < m_instruments.size(); ++i) {
            if (i > 0) {
                result.append(", ");
            }
            result.append(m_instruments[i]);
        }
    }
    else {
        result.append("all");
    }

    // Append Test Statistic
    if (m_has_ts) {
        result.append("\n"+gammalib::parformat("Test Statistic"));
        result.append(gammalib::str(ts()));
    }
    else if (m_tscalc) {
        result.append("\n"+gammalib::parformat("Test Statistic"));
        result.append("Computation requested");
    }

    // Append observation identifiers
    result.append("\n"+gammalib::parformat("Observation identifiers"));
    if (!m_ids.empty()) {
        for (int i = 0; i < m_ids.size(); ++i) {
            if (i > 0) {
                result.append(", ");
            }
            result.append(m_ids[i]);
        }
    }
    else {
        result.append("all");
    }

    // Append associations
    if (!m_associations.is_empty()) {
        result.append("\n"+gammalib::parformat("Associations"));
        for (int i = 0; i < m_associations.size(); ++i) {
            if (i > 0) {
                result.append(", ");
            }
            result.append(m_associations[i].name());
        }
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Read instrument scales from XML element
 *
 * @param[in] xml XML source element.
 *
 * Reads the instrument scale factors from a tag with the following format
 *
 *      <scaling>
 *         <instrument name="LAT" scale="1.0" min="0.1" max="10.0" value="1.0" free="0"/>
 *         <instrument name="CTA" scale="1.0" min="0.1" max="10.0" value="0.5" free="0"/>
 *      </scaling>
 *
 * The instrument name is case sensitive, but any leading and trailing white
 * space will be removed.
 *
 * If no scaling tag is found, all instrument scale factors will be cleared.
 ***************************************************************************/
void GModel::read_scales(const GXmlElement& xml)
{
    // Clear any existing scales
    m_scales.clear();

    // Continue only if there is a <scaling> tag
    if (xml.elements("scaling") != 0) {

        // Get pointer on first instrument scale factors
        const GXmlElement* scales = xml.element("scaling", 0);

        // Determine number of scale factors
        int nscales = scales->elements("instrument");

        // Read all scale factors
        for (int i = 0; i < nscales; ++i) {
            const GXmlElement* par = scales->element("instrument", i);
            GModelPar scale;
            scale.read(*par);
            scale.name(gammalib::strip_whitespace(par->attribute("name")));
            scale.has_grad(true);
            m_scales.push_back(scale);
            m_pars.push_back(&m_scales[m_scales.size()-1]);
        }

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write instrument scales into XML element
 *
 * @param[in] xml XML source element.
 *
 * @exception GException::invalid_value
 *            Invalid number of instrument tags found in XML element.
 *
 * If there are instrument scale factors then add a tag with the following
 * format to the XML element:
 *
 *      <scaling>
 *         <instrument name="LAT" scale="1" min="0.1" max="10" value="1.0" free="0"/>
 *         <instrument name="CTA" scale="1" min="0.1" max="10" value="0.5" free="0"/>
 *      </scaling>
 ***************************************************************************/
void GModel::write_scales(GXmlElement& xml) const
{
    // Continue only is scale factors are present
    if (!m_scales.empty()) {

        // Get number of instruments
        int num = m_scales.size();

        // Initialise scaling tag
        GXmlElement* scale = NULL;

        // If no <scaling> tag exists then add one now with the required
        // number of instruments ...
        if (xml.elements("scaling") == 0) {
            scale = xml.append("scaling");
            for (int i = 0; i < num; ++i) {
                scale->append(GXmlElement("instrument"));
            }
        }

        // ... otherwise get first tag
        else {
            scale = xml.element("scaling", 0);
        }

        // Verify that scaling tag has the required number of instruments
        if (scale->elements() != num) {
            std::string msg = "Number of "+gammalib::str(scale->elements())+
                              " scale elements in XML file does not correspond "
                              "to expected number of "+gammalib::str(num)+
                              " elements. Please verify the XML format.";
            throw GException::invalid_value(G_WRITE_SCALES, msg);
        }
        int npars = scale->elements("instrument");
        if (npars != num) {
            std::string msg = "Number of "+gammalib::str(npars)+" \"instrument\" "
                              "scale elements in XML file does not correspond to "
                              "expected number of "+gammalib::str(num)+
                              " elements. Please verify the XML format.";
            throw GException::invalid_value(G_WRITE_SCALES, msg);
        }

        // Write all instruments
        for (int i = 0; i < num; ++i) {

            // Get instrument element
            GXmlElement* inst = scale->element("instrument", i);

            // Set instrument name
            inst->attribute("name", m_scales[i].name());

            // Write instrument scaling factor
            m_scales[i].write(*inst);

        } // endfor: looped over all instruments

    }

    // Return
    return;
}
