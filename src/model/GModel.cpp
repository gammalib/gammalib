/***************************************************************************
 *              GModel.cpp - Abstract virtual model base class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModel.cpp
 * @brief GModel class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GModel.hpp"

/* __ Method name definitions ____________________________________________ */

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
 ***************************************************************************/
GModel::GModel(const GXmlElement& xml)
{
    // Initialise members
    init_members();

    // Set model name
    name(xml.attribute("name"));

    // Set instruments
    instruments(xml.attribute("instrument"));

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
 * @param[in] model Model which should be assigned.
 ***************************************************************************/
GModel& GModel::operator= (const GModel& model)
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


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

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
    for (int i = 0; i < inst.size(); ++i)
        m_instruments.push_back(toupper(strip_whitespace(inst[i])));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns instruments to which model applies
 *
 * @param[in] instruments String of instruments.
 *
 * Returns a comma separated list of instruments to which model applies. If
 * no instruments exist then an empty list is returned.
 ***************************************************************************/
std::string GModel::instruments(void) const
{
    // Initialise string
    std::string result;
    
    // Attach all instruments
    for (int i = 0; i < m_instruments.size(); ++i) {
        if (i > 0)
            result += ",";
        result += m_instruments[i];
    }

    // Return
    return result;
}


/***********************************************************************//**
 * @brief Verifies if model is valid for a given instrument
 *
 * @param[in] name Instrument name.
 *
 * Checks if specified instrument name is in list of applicable instruments.
 * If the list of applicable instruments is empty the model applies to all
 * possible instruments.
 ***************************************************************************/
bool GModel::isvalid(const std::string& name) const
{
    // Initialise validity
    bool valid = true;

    // Check if model applies to instrument
    if (m_instruments.size() > 0) {

        // Convert instrument name to upper case
        std::string uname = toupper(name);

        // Initialise validity flag
        valid = false;

        // Check if name is in list
        for (int i = 0; i < m_instruments.size(); ++i) {
            if (uname == m_instruments[i]) {
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

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model GModel members which should be copied.
 ***************************************************************************/
void GModel::copy_members(const GModel& model)
{
    // Copy attributes
    m_name        = model.m_name;
    m_instruments = model.m_instruments;

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
 * @brief Print model name and instrument
 ***************************************************************************/
std::string GModel::print_name_instrument(void) const
{
    // Initialise result string
    std::string result;

    // Append model
    result.append(parformat("Name")+name());
    result.append("\n"+parformat("Instruments"));
    if (m_instruments.size() > 0) {
        for (int i = 0; i < m_instruments.size(); ++i) {
            if (i > 0)
                result.append(", ");
            result.append(m_instruments[i]);
        }
    }
    else
        result.append("all");

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                                  Friends                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] model Model.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GModel& model)
{
     // Write spectrum in output stream
    os << model.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] model Model.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GModel& model)
{
    // Write spectrum into logger
    log << model.print();

    // Return logger
    return log;
}
