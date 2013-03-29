/***************************************************************************
 *              GModel.cpp - Abstract virtual model base class             *
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
#define G_ACCESS2                          "GModel::operator[](std::string&)"
#define G_AT                                               "GModel::at(int&)"
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
 * Construct model from XML element. The method extracts the "name",
 * "instrument" and "id" attributes from the model node.
 ***************************************************************************/
GModel::GModel(const GXmlElement& xml)
{
    // Initialise members
    init_members();

    // Set model name
    name(xml.attribute("name"));

    // Set instruments
    instruments(xml.attribute("instrument"));

    // Set observation identifiers
    ids(xml.attribute("id"));

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
 * @exception GException::par_not_found
 *            Parameter with specified name not found in container.
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
        throw GException::par_not_found(G_ACCESS2, name);
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
 * @exception GException::par_not_found
 *            Parameter with specified name not found in container.
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
        throw GException::par_not_found(G_ACCESS2, name);
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
 * @param[in] index Parameter index [0,...,size()-1].
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
        throw GException::out_of_range(G_AT, index, 0, size()-1);
    }

    // Return reference
    return *(m_pars[index]);
}


/***********************************************************************//**
 * @brief Returns reference to model parameter by index (const version)
 *
 * @param[in] index Parameter index [0,...,size()-1].
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
        throw GException::out_of_range(G_AT, index, 0, size()-1);
    }

    // Return reference
    return *(m_pars[index]);
}


/***********************************************************************//**
 * @brief Returns instruments to which model applies
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
 * list of strings. If the @p instrument string is empty the model is
 * considered to apply to all instruments.
 ***************************************************************************/
void GModel::instruments(const std::string& instruments)
{
    // Clear instruments vector
    m_instruments.clear();

    // Extract instruments
    std::vector<std::string> inst = split(instruments, ",");

    // Attach all instruments
    for (int i = 0; i < inst.size(); ++i) {
        m_instruments.push_back(toupper(strip_whitespace(inst[i])));
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns model scale factor for a given instrument
 *
 * @param[in] instrument Instrument.
 *
 * Returns the model scale factor for a given @p instrument. If the
 * @p instrument is not found, the method returns a scale factor of unity.
 ***************************************************************************/
GModelPar GModel::scale(const std::string& instrument) const
{
    // Convert instrument to upper case string
    std::string uinstrument = toupper(strip_whitespace(instrument));

    // Initialise unit scale factor
    GModelPar scale;
    scale.value(1.0);
    scale.name(uinstrument);
    scale.fix();

    // Search for instrument and recover scale factor if the instrument
    // has been found.
    for (int i = 0; i < m_scales.size(); ++i) {
        if (m_scales[i].name() == uinstrument) {
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
 * Sets the model parameter for a given instrument. The instrument name will
 * be determined from the model parameter name. The model parameter name
 * will be converted to upper case and any leading and trailing whitespace
 * will be stripped.
 *
 * If the instrument is not yet defined it will be appended to the list of
 * instruments.
 ***************************************************************************/
void GModel::scale(const GModelPar& par)
{
    // Convert instrument to upper case string
    std::string uinstrument = toupper(strip_whitespace(par.name()));

    // Search for instrument and copy the model parameter if the instrument
    // has been found. Make sure that the instrument name is in upper case.
    bool found = false;
    for (int i = 0; i < m_scales.size(); ++i) {
        if (m_scales[i].name() == uinstrument) {
            found       = true;
            m_scales[i] = par;
            m_scales[i].name(uinstrument);
            break;
        }
    }

    // If instrument has not been found then append it now to the list
    // of instruments.
    if (!found) {
        m_scales.push_back(par);
        m_scales[m_scales.size()-1].name(uinstrument);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns observation identifiers to which model applies
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
 * separated list of strings. If the observation identifier string is empty,
 * the model is considered to apply to all observation identifiers.
 ***************************************************************************/
void GModel::ids(const std::string& ids)
{
    // Clear observation identifier vector
    m_ids.clear();

    // Extract observation identifiers
    std::vector<std::string> id = split(ids, ",");

    // Attach all observation identifiers
    for (int i = 0; i < id.size(); ++i) {
        m_ids.push_back(toupper(strip_whitespace(id[i])));
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Verifies if model is valid for a given instrument and identifier
 *
 * @param[in] instrument Instrument name.
 * @param[in] id Observation identifier.
 * @return Validity flag
 *
 * Checks if specified instrument name and observation identifier is in list
 * of applicable instruments and identifiers. The check is case insensitive.
 *
 * If the list of applicable instruments is empty, the model applies to all
 * possible instruments. If the list of applicable observation identifiers
 * is empty, the model applies to all identifiers.
 ***************************************************************************/
bool GModel::isvalid(const std::string& instrument,
                     const std::string& id) const
{
    // Initialise validity
    bool valid = true;

    // Check if model applies to instrument
    if (!m_instruments.empty()) {

        // Convert instrument name to upper case
        std::string uinstrument = toupper(instrument);

        // Initialise validity flag
        valid = false;

        // Check if instrument is in list
        for (int i = 0; i < m_instruments.size(); ++i) {
            if (uinstrument == m_instruments[i]) {
                valid = true;
                break;
            }
        }

    }

    // Check if model applies to observation identifier
    if (valid && !m_ids.empty()) {

        // Convert observation identifier to upper case
        std::string uid = toupper(id);

        // Initialise validity flag
        valid = false;

        // Check if name is in list
        for (int i = 0; i < m_ids.size(); ++i) {
            if (uid == m_ids[i]) {
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
    m_name        = model.m_name;
    m_instruments = model.m_instruments;
    m_scales      = model.m_scales;
    m_ids         = model.m_ids;
    m_pars        = model.m_pars;

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
 * The  name of each instrument will be converted to upper case and any
 * leading or trailing whitespace will be removed. 
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
            scale.name(toupper(strip_whitespace(par->attribute("name"))));
            m_scales.push_back(scale);
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
 * @exception GException::model_invalid_parnum
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

        // Verify that scaling tag  has the required number of instruments
        if (scale->elements() != num || scale->elements("instrument") != num) {
            throw GException::model_invalid_parnum(G_WRITE_SCALES, *scale,
                  "Instrument scaling needs "+str(num)+" instrument tags.");
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
    result.append(parformat("Name")+name());

    // Append instruments
    result.append("\n"+parformat("Instruments"));
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

    // Append instrument scale factors
    result.append("\n"+parformat("Instrument scale factors"));
    if (!m_scales.empty()) {
        for (int i = 0; i < m_scales.size(); ++i) {
            if (i > 0) {
                result.append(", ");
            }
            result.append(m_scales[i].name());
            result.append("=");
            result.append(str(m_scales[i].value()));
        }
        result.append(", others unity");
    }
    else {
        result.append("unity");
    }

    // Append observation identifiers
    result.append("\n"+parformat("Observation identifiers"));
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

    // Return result
    return result;
}
