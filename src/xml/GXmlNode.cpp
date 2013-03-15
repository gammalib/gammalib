/***************************************************************************
 *                GXmlNode.cpp - Abstract XML node base class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GXmlNode.cpp
 * @brief Abstract XML node base class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GXmlNode.hpp"
#include "GXmlElement.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_APPEND                                "GXmlNode::append(GXmlNode*)"
#define G_CHILD                             "GXmlNode* GXmlNode::child(int&)"
#define G_ELEMENT1                        "GXmlNode* GXmlNode::element(int&)"
#define G_ELEMENT2          "GXmlNode* GXmlNode::element(std::string&, int&)"

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
GXmlNode::GXmlNode(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] node XML node.
 ***************************************************************************/
GXmlNode::GXmlNode(const GXmlNode& node)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(node);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GXmlNode::~GXmlNode(void)
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
 * @param[in] node XML node.
 * @return XML node.
 ***************************************************************************/
GXmlNode& GXmlNode::operator=(const GXmlNode& node)
{
    // Execute only if object is not identical
    if (this != &node) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(node);

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
 * @brief Append child node to node
 *
 * @param[in] node Child node.
 *
 * @exception GException::file_open_error
 *            Unable to open XML file (write access requested).
 ***************************************************************************/
void GXmlNode::append(GXmlNode* node)
{
    // Make sure that node is not a document (only the root document is
    // allowed to exist in a XML document
    if (node->type() == NT_DOCUMENT) {
        throw GException::xml_bad_node_type(G_APPEND, "GXmlDocument",
                          "Only the root node is of type GXmlDocument.");
    }

    // Append node
    m_nodes.push_back(node);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return number of children of XML node
 *
 * @return Number of children of XML node.
 *
 * Returns the number of child elements of the XML node. Child elements are
 * the elements that are directly contained by XML node.
 ***************************************************************************/
int GXmlNode::children(void) const
{
    // Return number of nodes
    return m_nodes.size();
}


/***********************************************************************//**
 * @brief Return pointer to child XML node
 *
 * @param[in] index Node index (0,...,children()-1).
 * @return Pointer to child of XML node.
 *
 * @exception GException::out_of_range
 *            Node index is out of range.
 *
 * Returns a pointer to the child number @p index of the XML node. An
 * exception will be thrown if the @index is not valid.
 ***************************************************************************/
GXmlNode* GXmlNode::child(const int& index) const
{
    // If index is outside boundary then throw an error
    if (index < 0 || index >= children()) {
        throw GException::out_of_range(G_CHILD, index, 0, children()-1);
    }

    // Return node
    return m_nodes[index];
}


/***********************************************************************//**
 * @brief Return number of GXMLElement children of node
 *
 * @return Number of GXMLElement children of XML node.
 *
 * Returns the number of GXMLElement child elements of the XML node.
 * GXMLElement child elements are nodes of type NT_ELEMENT.
 ***************************************************************************/
int GXmlNode::elements(void) const
{
    // Compute number of child elements in node
    int elements = 0;
    for (int i = 0; i < children(); ++i) {
        if (child(i)->type() == NT_ELEMENT) {
            elements++;
        }
    }

    // Return number of child elements
    return elements;
}


/***********************************************************************//**
 * @brief Return number of GXMLElement children with a given name
 *
 * @param[in] name Name of GXMLElement elements.
 * @return Number of GXMLElement children with a given @p name.
 *
 * Returns the number of GXMLElement child elements of the XML node that
 * have a given @p name. GXMLElement child elements are nodes of type
 * NT_ELEMENT.
 ***************************************************************************/
int GXmlNode::elements(const std::string& name) const
{
    // Compute number of child elements in node
    int elements = 0;
    for (int i = 0; i < children(); ++i) {
        if (child(i)->type() == NT_ELEMENT) {
            if (static_cast<GXmlElement*>(child(i))->name() == name) {
                elements++;
            }
        }
    }

    // Return number of child elements
    return elements;
}


/***********************************************************************//**
 * @brief Return pointer to GXMLElement child
 *
 * @param[in] index Node index (0,...,elements()-1).
 * @return Pointer to GXMLElement child.
 *
 * @exception GException::out_of_range
 *            Child element index is out of range.
 *
 * Returns a pointer to the child number @p index of the XML node. An
 * exception will be thrown if the @index is not valid.
 ***************************************************************************/
GXmlNode* GXmlNode::element(const int& index) const
{
    // If index is outside boundary then throw an error
    if (index < 0 || index >= elements()) {
        throw GException::out_of_range(G_ELEMENT1, index, 0, elements()-1);
    }

    // Compute number of child elements in node
    GXmlNode* element  = NULL;
    int       elements = 0;
    for (int i = 0; i < children(); ++i) {
        if (child(i)->type() == NT_ELEMENT) {
            if (elements == index) {
                element = child(i);
                break;
            }
            elements++;
        }
    }

    // Return child element
    return element;
}


/***********************************************************************//**
 * @brief Return pointer on GXMLElement child of a given name
 *
 * @param[in] name Name of child element.
 * @param[in] index Node index (0,...,elements()-1).
 *
 * @exception GException::xml_name_not_found
 *            Child element name not found.
 * @exception GException::out_of_range
 *            Child element index is out of range.
 *
 * Returns a pointer to the child number @p index with @p name of the XML
 * node. An exception will be thrown if the @index is not valid.
 ***************************************************************************/
GXmlNode* GXmlNode::element(const std::string& name, const int& index) const
{
    // Determine number of child elements
    int n = elements(name);

    // Signal if no children exist
    if (n < 1) {
        throw GException::xml_name_not_found(G_ELEMENT2, name);
    }

    // If index is outside boundary then throw an error
    if (index < 0 || index >= n) {
        throw GException::out_of_range(G_ELEMENT2, index, 0, n-1);
    }

    // Compute number of child elements in node
    GXmlNode* element  = NULL;
    int       elements = 0;
    for (int i = 0; i < children(); ++i) {
        if (child(i)->type() == NT_ELEMENT) {
            if (static_cast<GXmlElement*>(child(i))->name() == name) {
                if (elements == index) {
                    element = child(i);
                    break;
                }
                elements++;
            }
        }
    }

    // Return child element
    return element;
}


/***********************************************************************//**
 * @brief Print XML node in string
 *
 * @return String containing XML information.
 ***************************************************************************/
std::string GXmlNode::print(void) const
{
    // Set result string
    std::string result = print(0);

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
void GXmlNode::init_members(void)
{
    // Initialise members
    m_nodes.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] node XML node.
 ***************************************************************************/
void GXmlNode::copy_members(const GXmlNode& node)
{
    // Copy nodes
    m_nodes.clear();
    for (int i = 0; i < node.m_nodes.size(); ++i) {
        m_nodes.push_back((node.m_nodes[i]->clone()));
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 *
 * As container classes that hold pointers need to handle themselves the
 * proper deallocation of memory, we loop here over all pointers and make
 * sure that we deallocate the associated nodes.
 ***************************************************************************/
void GXmlNode::free_members(void)
{
    // Free memory for all nodes
    for (int i = 0; i < m_nodes.size(); ++i) {
        delete m_nodes[i];
        m_nodes[i] = NULL;
    }

    // Return
    return;
}
