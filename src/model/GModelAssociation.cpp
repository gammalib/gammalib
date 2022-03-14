/***************************************************************************
 *              GModelAssociation.cpp - Model association class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020-2022 by Juergen Knoedlseder                         *
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
 * @file GModelAssociation.cpp
 * @brief Model association class implementation
 * @author Juergen Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GXmlElement.hpp"
#include "GModelAssociation.hpp"

/* __ Constants __________________________________________________________ */

/* __ Method name definitions ____________________________________________ */
#define G_VALUE                      "GModelAssociation::value(std::string&)"
#define G_ERROR                      "GModelAssociation::error(std::string&)"
#define G_PROPERTY  "GModelAssociation::property(std::string&, std::string&,"\
                                                             " std::string&)"
#define G_GET_PROPERTY_XML             "GModelAssociation::get_property_xml("\
                                                "GXmlElement&, std::string&)"

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
GModelAssociation::GModelAssociation(void)
{
    // Initialise private members
    init_members();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Association name constructor
 *
 * @param[in] name Association name.
 *
 * Construct a model association from an association name.
 ***************************************************************************/
GModelAssociation::GModelAssociation(const std::string& name)
{
    // Initialise private members
    init_members();

    // Set association name
    this->name(name);

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML element constructor
 *
 * @param[in] xml XML element.
 *
 * Construct a model association from an XML element.
 ***************************************************************************/
GModelAssociation::GModelAssociation(const GXmlElement& xml)
{
    // Initialise private members
    init_members();

    // Read model association from XML element
    read(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] association Model association.
 ***************************************************************************/
GModelAssociation::GModelAssociation(const GModelAssociation& association)
{ 
    // Initialise private members
    init_members();

    // Copy members
    copy_members(association);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GModelAssociation::~GModelAssociation(void)
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
 * @param[in] association Model association.
 * @return Model association.
 ***************************************************************************/
GModelAssociation& GModelAssociation::operator=(const GModelAssociation& association)
{ 
    // Execute only if object is not identical
    if (this != &association) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(association);

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
 * @brief Clear model association
 ***************************************************************************/
void GModelAssociation::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();
    
    // Return
    return; 
}


/***********************************************************************//**
 * @brief Clone model association
 *
 * @return Pointer to deep copy of model association.
 ***************************************************************************/
GModelAssociation* GModelAssociation::clone(void) const
{
    // Clone this model association
    return new GModelAssociation(*this);
}


/***********************************************************************//**
 * @brief Return property value
 *
 * @param[in] name Property name
 * @return Property value.
 *
 * @exception GException::invalid_argument
 *            Property name not found
 *
 * Returns the value of the property with the specified @p name.
 ***************************************************************************/
const std::string& GModelAssociation::value(const std::string& name) const
{
    // Get model index
    int index = get_index(name);

    // Throw exception if property name was not found
    if (index == -1) {
        std::string msg = "Property \""+name+"\" not found in association.";
        throw GException::invalid_argument(G_VALUE, msg);
    }

    // Return property value
    return (m_values[index]);
}


/***********************************************************************//**
 * @brief Return property error
 *
 * @param[in] name Property name
 * @return Property error.
 *
 * @exception GException::invalid_argument
 *            Property name not found
 *
 * Returns the error of the property with the specified @p name.
 ***************************************************************************/
const std::string& GModelAssociation::error(const std::string& name) const
{
    // Get model index
    int index = get_index(name);

    // Throw exception if property name was not found
    if (index == -1) {
        std::string msg = "Association property \""+name+"\" not found.";
        throw GException::invalid_argument(G_ERROR, msg);
    }

    // Return property error
    return (m_errors[index]);
}


/***********************************************************************//**
* @brief Set property value and (optionally) error
*
* @param[in] name Property name
* @param[in] value Property value
* @param[in] error Property error
*
* @exception GException::invalid_argument
*            Property with specified @p name exists aleady
*
* Sets the value and optionally the error of the property with the
* specified @p name.
***************************************************************************/
void GModelAssociation::property(const std::string& name,
                                 const std::string& value,
                                 const std::string& error)
{
    // Throw an exception if property exists already
    if (get_index(name) != -1) {
        std::string msg = "Association property \""+name+"\" exists already.";
        throw GException::invalid_argument(G_PROPERTY, msg);
    }

    // Push name, value and error in property lists
    m_names.push_back(name);
    m_values.push_back(value);
    m_errors.push_back(error);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read model association from XML document
 *
 * @param[in] xml XML element.
 *
 * Read model association from the XML element. The XML element is expected
 * to have the following structure
 *
 *     <association name="Crab">
 *         <property name="RA" value="83.6331"/>
 *         <property name="DEC" value="22.0145"/>
 *         <property name="distance" value="0.0123"/>
 *         <property name="probability" value="0.978"/>
 *     </association>
 *
 * Properties need to have unique names. If properties with identical names
 * are encountered, and exception is thrown.
 ***************************************************************************/
void GModelAssociation::read(const GXmlElement& xml)
{
    // Initialise properties
    m_names.clear();
    m_values.clear();
    m_errors.clear();

    // Set association name
    name(xml.attribute("name"));

    // Get number of properties
    int n = xml.elements("property");

    // Loop over all properties
    for (int i = 0; i < n; ++i) {

        // Get pointer on property
        const GXmlElement* property = xml.element("property", i);

        // Read name, value and error. If one of the attributes does not
        // exist an empty string will be returned
        std::string name  = property->attribute("name");
        std::string value = property->attribute("value");
        std::string error = property->attribute("error");

        // Set property
        this->property(name, value, error);

    } // endfor: looped over properties

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model association into XML element
 *
 * @param[in] xml XML element.
 *
 * Write model association into the XML element. The following structure
 * will be written
 *
 *     <association name="Crab">
 *         <property name="RA" value="83.6331"/>
 *         <property name="DEC" value="22.0145"/>
 *         <property name="distance" value="0.0123"/>
 *         <property name="probability" value="0.978"/>
 *     </association>
 *
 ***************************************************************************/
void GModelAssociation::write(GXmlElement& xml) const
{
    // Write association name
    xml.attribute("name", name());

    // Loop over all association properties
    for (int i = 0; i < size(); ++i) {
        GXmlElement* property = get_property_xml(xml, m_names[i]);
        property->attribute("name", m_names[i]);
        property->attribute("value", m_values[i]);
        if (!m_errors[i].empty()) {
            property->attribute("error", m_errors[i]);
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print model association
 *
 * @param[in] chatter Chattiness.
 * @return String containing model association information.
 ***************************************************************************/
std::string GModelAssociation::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelAssociation ===");

        // Append association name
        result.append("\n"+gammalib::parformat("Name")+name());

        // Append association properties
        result.append("\n"+gammalib::parformat("Number of properties"));
        result.append(gammalib::str(size()));
        for (int i = 0; i < size(); ++i) {
            result.append("\n"+gammalib::parformat(m_names[i])+m_values[i]);
            if (m_errors[i] != "") {
                result.append(" +/- "+m_errors[i]);
            }
        }

    } // endif: chatter was not silent

    // Return
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModelAssociation::init_members(void)
{
    // Initialise members
    m_name.clear();
    m_names.clear();
    m_values.clear();
    m_errors.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] association Model association.
 ***************************************************************************/
void GModelAssociation::copy_members(const GModelAssociation& association)
{
    // Copy members
    m_name   = association.m_name;
    m_names  = association.m_names;
    m_values = association.m_values;
    m_errors = association.m_errors;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelAssociation::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Return property index by name
 *
 * @param[in] name Property name.
 * @return Property index (-1 if not found)
 *
 * Returns the index of the property with the specified @p name. If no
 * property with the specified @p name is found the method returns -1.
 ***************************************************************************/
int GModelAssociation::get_index(const std::string& name) const
{
    // Initialise index
    int index = -1;

    // Search property names for the specified name
    for (int i = 0; i < size(); ++i) {
        if (m_names[i] == name) {
            index = i;
            break;
        }
    }

    // Return index
    return index;
}


/***********************************************************************//**
 * @brief Return pointer to property with given name in XML element
 *
 * @param[in] xml XML element.
 * @param[in] name Property name.
 * @return Pointer to property XML element.
 *
 * @exception GException::invalid_value
 *            Invalid XML format encountered.
 *
 * Returns pointer to property with given @p name in XML element. If the
 * @p name is not found, a property with the given @p name is added.
 *
 * The function checks for multiple occurences of a property and throws an
 * exception in case that more than one property with a given name is found.
 ***************************************************************************/
GXmlElement* GModelAssociation::get_property_xml(GXmlElement&       xml,
                                                 const std::string& name) const
{
    // Initialize XML element pointer
    GXmlElement* property = NULL;

    // Number of elements
    int number = 0;

    // Get number of elements in XML element
    int n = xml.elements("property");

    // Search for property with given name
    for (int i = 0; i < n; ++i) {
        GXmlElement* element = xml.element("property", i);
        if (element->attribute("name") == name) {
            property = element;
            number++;
        }
    }

    // Create property if none was found
    if (number == 0) {
        property = static_cast<GXmlElement*>(xml.append(GXmlElement("property name=\""+name+"\"")));
        number++;
    }

    // Throw case dependent exception
    if (number < 1) {
        std::string msg = "Property \""+name+"\" not found in XML element."
                          " Please verify the XML format.";
        throw GException::invalid_value(G_GET_PROPERTY_XML, msg);
    }
    else if (number > 1) {
        std::string msg = "Property \""+name+"\" found "+
                          gammalib::str(number)+" times in XML element."
                          " Please verify the XML format.";
        throw GException::invalid_value(G_GET_PROPERTY_XML, msg);
    }

    // Return
    return property;
}
