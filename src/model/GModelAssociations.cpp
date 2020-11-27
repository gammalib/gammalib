/***************************************************************************
 *       GModelAssociations.cpp - Model associations container class       *
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
 * @file GModelAssociations.cpp
 * @brief Model associations container class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GXmlElement.hpp"
#include "GModelAssociations.hpp"
#include "GModelAssociation.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ACCESS               "GModelAssociations::operator[](std::string&)"
#define G_AT                                   "GModelAssociations::at(int&)"
#define G_SET1            "GModelAssociations::set(int&, GModelAssociation&)"
#define G_SET2    "GModelAssociations::set(std::string&, GModelAssociation&)"
#define G_APPEND             "GModelAssociations::append(GModelAssociation&)"
#define G_INSERT1      "GModelAssociations::insert(int&, GModelAssociation&)"
#define G_INSERT2                 "GModelAssociations::insert(std::string&, "\
                                                        "GModelAssociation&)"
#define G_REMOVE1                          "GModelAssociations::remove(int&)"
#define G_REMOVE2                  "GModelAssociations::remove(std::string&)"
#define G_EXTEND            "GModelAssociations::extend(GModelAssociations&)"
#define G_READ                              "GModelAssociations::read(GXml&)"
#define G_GET_ASSOCIATION_XML      "GModelAssociations::get_association_xml("\
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
GModelAssociations::GModelAssociations(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML element constructor
 *
 * @param[in] xml XML element.
 *
 * Construct a model associations container from an XML element.
 ***************************************************************************/
GModelAssociations::GModelAssociations(const GXmlElement& xml)
{
    // Initialise private members
    init_members();

    // Read model associations container from XML element
    read(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] associations Model associations container.
 ***************************************************************************/
GModelAssociations::GModelAssociations(const GModelAssociations& associations)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(associations);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GModelAssociations::~GModelAssociations(void)
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
 * @param[in] associations Model associations container.
 * @return Model associations container.
 ***************************************************************************/
GModelAssociations& GModelAssociations::operator=(const GModelAssociations& associations)
{
    // Execute only if object is not identical
    if (this != &associations) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(associations);

    } // endif: object was not identical

    // Return
    return *this;
}


/***********************************************************************//**
 * @brief Return reference to model association
 *
 * @param[in] name Model association name.
 *
 * @exception GException::invalid_argument
 *            Model association with specified @p name not found in container.
 *
 * Returns a reference to the model with the specified @p name.
 ***************************************************************************/
GModelAssociation& GModelAssociations::operator[](const std::string& name)
{
    // Return association using const method
    return const_cast<GModelAssociation&>((*static_cast<const GModelAssociations*>(this))[name]);
}


/***********************************************************************//**
 * @brief Return reference to model association (const version)
 *
 * @param[in] name Model association name.
 *
 * @exception GException::invalid_argument
 *            Model association with specified @p name not found in container.
 *
 * Returns a const reference to the model with the specified @p name.
 ***************************************************************************/
const GModelAssociation& GModelAssociations::operator[](const std::string& name) const
{
    // Get model index
    int index = get_index(name);

    // Throw exception if model name was not found
    if (index == -1) {
        std::string msg = "Model association \""+name+"\" not found in "
                          "container.";
        throw GException::invalid_argument(G_ACCESS, name);
    }

    // Return pointer
    return m_associations[index];
}


/*==========================================================================
 =                                                                         =
 =                            Public methods                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear object
 *
 * Removes all model associations from the container.
 ***************************************************************************/
void GModelAssociations::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 *
 * @return Pointer to deep copy of model associations container
 *
 * Makes a deep copy of the model associations container instance.
 ***************************************************************************/
GModelAssociations* GModelAssociations::clone(void) const
{
    return new GModelAssociations(*this);
}


/***********************************************************************//**
 * @brief Return reference to model association
 *
 * @param[in] index Model association index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Model association index is out of range.
 *
 * Returns a reference to the model association with the specified @p index.
 ***************************************************************************/
GModelAssociation& GModelAssociations::at(const int& index)
{
    // Return association using const method
    return const_cast<GModelAssociation&>(static_cast<const GModelAssociations*>(this)->at(index));
}


/***********************************************************************//**
 * @brief Return reference to model association (const version)
 *
 * @param[in] index Model association index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Model association index is out of range.
 *
 * Returns a const reference to the model association with the specified
 * @p index.
 ***************************************************************************/
const GModelAssociation& GModelAssociations::at(const int& index) const
{
    // Compile option: raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Model association index", index, size());
    }
    #endif

    // Return pointer
    return m_associations[index];
}


/***********************************************************************//**
 * @brief Append model association to container
 *
 * @param[in] association Model association.
 * @return Reference to appended model association.
 *
 * @exception GException::invalid_value
 *            Name of model association exists already in container.
 *
 * Appends model association to the container.
 ***************************************************************************/
GModelAssociation& GModelAssociations::append(const GModelAssociation& association)
{
    // Check if a model association with specified name does not yet exist
    int inx = get_index(association.name());
    if (inx != -1) {
        std::string msg =
            "Attempt to append model association \""+association.name()+"\" "
            "to association container, but an association with the same name "
            "exists already at index "+gammalib::str(inx)+" in the container. "
            "Every model association in the container needs a unique name.";
        throw GException::invalid_value(G_APPEND, msg);
    }

    // Append model association
    m_associations.push_back(association);

    // Return reference to model association
    return (m_associations[size()-1]);
}


/***********************************************************************//**
 * @brief Insert model association into container
 *
 * @param[in] index Model association index [0,...,size()-1].
 * @param[in] association Model association.
 * @return Reference to inserted model association.
 *
 * @exception GException::out_of_range
 *            Model association index is out of range.
 * @exception GException::invalid_value
 *            Name of model association exists already in container.
 *
 * Inserts a @p model association into the container before the model
 * association with the specified @p index.
 ***************************************************************************/
GModelAssociation& GModelAssociations::insert(const int&               index,
                                              const GModelAssociation& association)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (is_empty()) {
        if (index > 0) {
            throw GException::out_of_range(G_INSERT1, "Model association index",
                                           index, size());
        }
    }
    else {
        if (index < 0 || index >= size()) {
            throw GException::out_of_range(G_INSERT1, "Model association index",
                                           index, size());
        }
    }
    #endif

    // Check if a model association with specified name does not yet exist
    int inx = get_index(association.name());
    if (inx != -1) {
        std::string msg =
            "Attempt to insert model association \""+association.name()+"\" "
            "in container before index "+gammalib::str(index)+", but an "
            "association with the same name exists already at index "+
            gammalib::str(inx)+" in the container. Every model association "
            "in the container needs a unique name.";
        throw GException::invalid_value(G_INSERT1, msg);
    }

    // Inserts deep copy of model association
    m_associations.insert(m_associations.begin()+index, association);

    // Return reference to model association
    return (m_associations[index]);
}


/***********************************************************************//**
 * @brief Insert model association into container
 *
 * @param[in] name Model association name.
 * @param[in] association Model association.
 * @return Reference to inserted model association.
 *
 * @exception GException::invalid_argument
 *            Model association with specified name not found in container.
 * @exception GException::invalid_value
 *            Name of model association exists already in container.
 *
 * Inserts a @p model association into the container before the model
 * association with the specified @p name.
 ***************************************************************************/
GModelAssociation& GModelAssociations::insert(const std::string&       name,
                                              const GModelAssociation& association)
{
    // Get association index
    int index = get_index(name);

    // Throw exception if parameter name was not found
    if (index == -1) {
        std::string msg = "Model association \""+name+"\" not found in "
                          "container.";
        throw GException::invalid_argument(G_INSERT2, name);
    }

    // Check if a model with specified name does not yet exist
    int inx = get_index(association.name());
    if (inx != -1) {
        std::string msg =
            "Attempt to insert model association \""+association.name()+"\" "
            "in container before index "+gammalib::str(index)+", but an "
            "association with the same name exists already at index "+
            gammalib::str(inx)+" in the container. Every model association "
            "in the container needs a unique name.";
        throw GException::invalid_value(G_INSERT2, msg);
    }

    // Inserts deep copy of model association
    m_associations.insert(m_associations.begin()+index, association);

    // Return reference to model association
    return (m_associations[index]);
}


/***********************************************************************//**
 * @brief Remove model association from container
 *
 * @param[in] index Model association index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Model association index is out of range.
 *
 * Remove model association of specified @p index from container.
 ***************************************************************************/
void GModelAssociations::remove(const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_REMOVE1, "Model association index",
                                       index, size());
    }
    #endif

    // Erase model association component from container
    m_associations.erase(m_associations.begin() + index);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Remove model association from container
 *
 * @param[in] name Model association name.
 *
 * @exception GException::invalid_argument
 *            Model association with specified name not found in container.
 *
 * Remove model association of specified @p name from container.
 ***************************************************************************/
void GModelAssociations::remove(const std::string& name)
{
    // Get parameter index
    int index = get_index(name);

    // Throw exception if parameter name was not found
    if (index == -1) {
        std::string msg = "Model association \""+name+"\" not found in "
                          "container.";
        throw GException::invalid_argument(G_REMOVE2, name);
    }

    // Erase model association component from container
    m_associations.erase(m_associations.begin() + index);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append model association container
 *
 * @param[in] associations Model association container.
 *
 * Append model association container to the container.
 ***************************************************************************/
void GModelAssociations::extend(const GModelAssociations& associations)
{
    // Do nothing if model association container is empty
    if (!associations.is_empty()) {

        // Get size. Note that we extract the size first to avoid an
        // endless loop that arises when a container is appended to
        // itself.
        int num = associations.size();

        // Reserve enough space
        reserve(size() + num);

        // Loop over all model association components and append pointers
        // to deep copies
        for (int i = 0; i < num; ++i) {

            // Check if model association name does not yet exist
            int inx = get_index(associations[i].name());
            if (inx != -1) {
                std::string msg =
                    "Attempt to append model association \""+
                    associations[i].name()+"\" to container, but an "
                    "association with the same name exists already at "
                    "index "+gammalib::str(inx)+" in the container. Every "
                    "model association in the container needs a unique name.";
                throw GException::invalid_value(G_EXTEND, msg);
            }

            // Append model association to container
            m_associations.push_back(associations[i]);

        } // endfor: looped over all model associations

    } // endif: model association container was not empty
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Signals if model association name exists
 *
 * @param[in] name Model association name.
 * @return True if model association with specified @p name exists.
 *
 * Searches all model association names for a match with the specified
 * @p name. If the specified name has been found, true is returned.
 ***************************************************************************/
bool GModelAssociations::contains(const std::string& name) const
{
    // Get model association index
    int index = get_index(name);

    // Return
    return (index != -1);
}


/***********************************************************************//**
 * @brief Read model associations from XML document
 *
 * @param[in] xml XML element.
 *
 * Read model associations from the XML element. The XML element is expected
 * to have the following structure
 *
 *     <associations>
 *         <association name="Crab">
 *             <property name="RA" value="83.6331"/>
 *             <property name="DEC" value="22.0145"/>
 *             <property name="distance" value="0.0123"/>
 *             <property name="probability" value="0.978"/>
 *         </association>
 *     </associations>
 *
 * This method does nothing if no "associations" tag is present in the XML
 * element.
 ***************************************************************************/
void GModelAssociations::read(const GXmlElement& xml)
{
    // Initialise associations
    m_associations.clear();

    // Continue only if there is an <associations> tag
    if (xml.elements("associations") > 0) {

        // Get pointer on associations
        const GXmlElement* associations = xml.element("associations", 0);

        // Determine number of associations
        int n = associations->elements("association");

        // Loop over all associations
        for (int i = 0; i < n; ++i) {

            // Get pointer on association
            const GXmlElement* element = associations->element("association", i);

            // Construct association from XML element
            GModelAssociation association(*element);

            // Append assocition to container
            append(association);

        } // endfor: looped over all associations

    } // endif: an association tag was present

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write models into XML element
 *
 * @param[in] xml XML element.
 *
 * Write model associations into the first "associations" element that is
 * found in the XML element. In case that no "associations" element exists,
 * one is added to the XML element. The written XML element has the
 * following structure
 *
 *     <associations>
 *         <association name="Crab">
 *             <property name="RA" value="83.6331"/>
 *             <property name="DEC" value="22.0145"/>
 *             <property name="distance" value="0.0123"/>
 *             <property name="probability" value="0.978"/>
 *         </association>
 *     </associations>
 *
 * This method does nothing if there are no associations in the instance.
 ***************************************************************************/
void GModelAssociations::write(GXmlElement& xml) const
{
    // Continue only if there are associations
    if (!is_empty()) {

        // If there is no associations tag then append one
        if (xml.elements("associations") == 0) {
            xml.append(GXmlElement("associations"));
        }

        // Get pointer on model associations
        GXmlElement* associations = xml.element("associations", 0);

        // Write all model associations into the XML element
        for (int i = 0; i < size(); ++i) {
            GXmlElement* association = get_association_xml(*associations,
                                       m_associations[i].name());
            m_associations[i].write(*association);
        }

    } // endif: there were associations to write

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print model associations
 *
 * @param[in] chatter Chattiness.
 * @return String containing model association container information.
 *
 * Prints all model associations into a string.
 ***************************************************************************/
std::string GModelAssociations::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelAssociations ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of associations"));
        result.append(gammalib::str(size()));
 
        // Append associations
        for (int i = 0; i < size(); ++i) {
            result.append("\n"+m_associations[i].print(chatter));
        }

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
void GModelAssociations::init_members(void)
{
    // Initialise members
    m_associations.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] associations Model association container.
 ***************************************************************************/
void GModelAssociations::copy_members(const GModelAssociations& associations)
{
    // Copy members
    m_associations = associations.m_associations;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelAssociations::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Return model association index by name
 *
 * @param[in] name Model ssociation name.
 * @return Model association index (-1 if not found)
 *
 * Returns model association index based on the specified @p name. If no
 * model association with the specified @p name is found the method returns
 * -1.
 ***************************************************************************/
int GModelAssociations::get_index(const std::string& name) const
{
    // Initialise index
    int index = -1;

    // Search model association with specified name
    for (int i = 0; i < size(); ++i) {
        if (m_associations[i].name() == name) {
            index = i;
            break;
        }
    }

    // Return index
    return index;
}


/***********************************************************************//**
 * @brief Return pointer to model association with given name in XML element
 *
 * @param[in] xml XML element.
 * @param[in] name Model association name.
 * @return Pointer to model association XML element.
 *
 * @exception GException::invalid_value
 *            Invalid XML format encountered.
 *
 * Returns pointer to model association with given @p name in XML element.
 * If the @p name is not found, a model association with the given @p name
 * is added.
 *
 * The function checks for multiple occurences of a model association and
 * throws an exception in case that more than one association with a given
 * name is found.
 ***************************************************************************/
GXmlElement* GModelAssociations::get_association_xml(GXmlElement&       xml,
                                                     const std::string& name) const
{
    // Initialize XML element pointer
    GXmlElement* association = NULL;

    // Number of elements
    int number = 0;

    // Get number of elements in XML element
    int n = xml.elements("association");

    // Search for property with given name
    for (int i = 0; i < n; ++i) {
        GXmlElement* element = xml.element("association", i);
        if (element->attribute("name") == name) {
            association = element;
            number++;
        }
    }

    // Create property if none was found
    if (number == 0) {
        association = static_cast<GXmlElement*>(xml.append(GXmlElement("association name=\""+name+"\"")));
        number++;
    }

    // Throw case dependent exception
    if (number < 1) {
        std::string msg = "Model association \""+name+"\" not found in XML "
                          " element. Please verify the XML format.";
        throw GException::invalid_value(G_GET_ASSOCIATION_XML, msg);
    }
    else if (number > 1) {
        std::string msg = "Model association \""+name+"\" found "+
                          gammalib::str(number)+" times in XML element."
                          " Please verify the XML format.";
        throw GException::invalid_value(G_GET_ASSOCIATION_XML, msg);
    }

    // Return
    return association;
}
