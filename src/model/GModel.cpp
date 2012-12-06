/***************************************************************************
 *              GModel.cpp - Abstract virtual model base class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2012 by Juergen Knoedlseder                         *
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
#define G_ACCESS1                                  "GModel::operator[](int&)"
#define G_ACCESS2                          "GModel::operator[](std::string&)"

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
 * @brief Returns reference to model parameter by index
 *
 * @param[in] index Parameter index [0,...,size()-1].
 * @return Reference to model parameter.
 *
 * @exception GException::out_of_range
 *            Parameter index is out of range.
 ***************************************************************************/
GModelPar& GModel::operator[](const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_ACCESS1, index, 0, size()-1);
    }
    #endif

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
 ***************************************************************************/
const GModelPar& GModel::operator[](const int& index) const
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_ACCESS1, index, 0, size()-1);
    }
    #endif

    // Return reference
    return *(m_pars[index]);
}


/***********************************************************************//**
 * @brief Returns reference to model parameter by name
 *
 * @param[in] name Parameter name.
 * @return Reference to model parameter.
 *
 * @exception GException::par_not_found
 *            Parameter with specified name not found in container.
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
 * list of strings. If the instrument string is empty the model is considered
 * to apply to all instruments.
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
 * @brief Print model attributes
 ***************************************************************************/
std::string GModel::print_attributes(void) const
{
    // Initialise result string
    std::string result;

    // Append model
    result.append(parformat("Name")+name());
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
