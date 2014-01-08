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
#define G_ACCESS                                 "GXmlNode::operator[](int&)"
#define G_SET                                "GXmlNode::set(int&, GXmlNode&)"
#define G_APPEND                                "GXmlNode::append(GXmlNode*)"
#define G_INSERT                          "GXmlNode::insert(int&, GXmlNode&)"
#define G_REMOVE                                     "GXmlNode::remove(int&)"
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


/***********************************************************************//**
 * @brief Return pointer to XML child node
 *
 * @param[in] index Child node index [0,...,size()-1].
 * @return Pointer to XML child node at @p index.
 *
 * @exception GException::out_of_range
 *            Child node index is out of range.
 *
 * Returns a pointer to the XML child node with the specified @p index.
 ***************************************************************************/
GXmlNode* GXmlNode::operator[](const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_ACCESS, index, 0, size()-1);
    }
    #endif

    // Return pointer
    return m_nodes[index];
}


/***********************************************************************//**
 * @brief Return pointer to XML child node (const version)
 *
 * @param[in] index Child node index [0,...,size()-1].
 * @return Pointer to XML child node at @p index.
 *
 * @exception GException::out_of_range
 *            Child node index is out of range.
 *
 * Returns a const pointer to the XML child node with the specified @p index.
 ***************************************************************************/
const GXmlNode* GXmlNode::operator[](const int& index) const
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_ACCESS, index, 0, size()-1);
    }
    #endif

    // Return pointer
    return m_nodes[index];
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Set XML child node
 *
 * @param[in] index Child node index [0,...,size()-1].
 * @param[in] node XML child node.
 * @return Pointer to deep copy of child node
 *
 * @exception GException::xml_bad_node_type
 *            Not allowed to append document node.
 * @exception GException::out_of_range
 *            Child node index is out of range.
 *
 * Set XML child node. A deep copy of the node will be made and the pointer
 * to this node will be stored.
 ***************************************************************************/
GXmlNode* GXmlNode::set(const int& index, const GXmlNode& node)
{
    // Make sure that node is not a document (only the root document is
    // allowed to exist in a XML document
    if (node.type() == NT_DOCUMENT) {
        throw GException::xml_bad_node_type(G_SET, "GXmlDocument",
                          "Only the root node is of type GXmlDocument.");
    }

    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_SET, index, 0, size()-1);
    }
    #endif

    // Delete any existing node
    if (m_nodes[index] != NULL) delete m_nodes[index];

    // Assign new child node by cloning
    m_nodes[index] = node.clone();

    // Return pointer
    return m_nodes[index];
}


/***********************************************************************//**
 * @brief Append XML child node
 *
 * @param[in] node XML child node.
 * @return Pointer to appended child node
 *
 * @exception GException::xml_bad_node_type
 *            Not allowed to append document node.
 *
 * Appends XML child node by making a deep copy of the node and storing its
 * pointer.
 ***************************************************************************/
GXmlNode* GXmlNode::append(const GXmlNode& node)
{
    // Make sure that node is not a document (only the root document is
    // allowed to exist in a XML document
    if (node.type() == NT_DOCUMENT) {
        throw GException::xml_bad_node_type(G_APPEND, "GXmlDocument",
                          "Only the root node is of type GXmlDocument.");
    }

    // Clone child node
    GXmlNode* ptr = node.clone();

    // Append deep copy of child node
    m_nodes.push_back(ptr);

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Append XML element child node
 *
 * @param[in] segment XML child node.
 * @return Pointer to appended child node
 *
 * Appends XML element that is constructed from a text @p segment. The text
 * segment is parsed and the element name and attributes are extracted using
 * the GXmlElement::parse_start() method. The method returns a pointer to the
 * XML element child node that has been appended.
 ***************************************************************************/
GXmlElement* GXmlNode::append(const std::string& segment)
{
    // Create a new XML child element from the text segment
    GXmlElement* ptr = new GXmlElement(segment);

    // Append child element to container
    m_nodes.push_back(ptr);

    // Return child element
    return ptr;
}


/***********************************************************************//**
 * @brief Insert XML child node
 *
 * @param[in] index Child node index [0,...,size()-1].
 * @param[in] node XML child node.
 * @return Pointer to inserted child node
 *
 * @exception GException::xml_bad_node_type
 *            Not allowed to append document node.
 * @exception GException::out_of_range
 *            Child node index is out of range.
 *
 * Inserts the XML child @p node before the node with the specified @p index.
 * A deep copy of the node will be made and the pointer to this node will be
 * stored.
 ***************************************************************************/
GXmlNode* GXmlNode::insert(const int& index, const GXmlNode& node)
{
    // Make sure that node is not a document (only the root document is
    // allowed to exist in a XML document
    if (node.type() == NT_DOCUMENT) {
        throw GException::xml_bad_node_type(G_INSERT, "GXmlDocument",
                          "Only the root node is of type GXmlDocument.");
    }

    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (is_empty()) {
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

    // Clone child node
    GXmlNode* ptr = node.clone();

    // Inserts deep copy of child node
    m_nodes.insert(m_nodes.begin()+index, ptr);

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Remove XML child node
 *
 * @param[in] index Child node index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Child node index is out of range.
 *
 * Remove XML child node at @p index.
 ***************************************************************************/
void GXmlNode::remove(const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_REMOVE, index, 0, size()-1);
    }
    #endif

    // Delete node
    delete m_nodes[index];
    
    // Erase child node from container
    m_nodes.erase(m_nodes.begin() + index);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Append all XML child nodes from another XML node
 *
 * @param[in] node XML node.
 *
 * Append all XML child nodes found in @p node to the actual object. Nodes
 * are copied deeply so that they live now on their on in the actual object.
 ***************************************************************************/
void GXmlNode::extend(const GXmlNode& node)
{
    // Do nothing if node container is empty
    if (!node.is_empty()) {

        // Get size. Note that we extract the size first to avoid an
        // endless loop that arises when a container is appended to
        // itself.
        int num = node.size();

        // Reserve enough space
        reserve(size() + num);

        // Loop over all child nodes and append pointers to deep copies 
        for (int i = 0; i < num; ++i) {
            m_nodes.push_back(node[i]->clone());
        }

    } // endif: node container was not empty
    
    // Return
    return;
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
    for (int i = 0; i < m_nodes.size(); ++i) {
        if (m_nodes[i]->type() == NT_ELEMENT) {
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
    for (int i = 0; i < m_nodes.size(); ++i) {
        if (m_nodes[i]->type() == NT_ELEMENT) {
            if (static_cast<GXmlElement*>(m_nodes[i])->name() == name) {
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
 * @param[in] index Node index [0,...,elements()-1].
 * @return Pointer to GXMLElement child.
 *
 * @exception GException::out_of_range
 *            Child element index is out of range.
 *
 * Returns a pointer to the child number @p index of the XML node. An
 * exception will be thrown if the @p index is not valid.
 ***************************************************************************/
GXmlElement* GXmlNode::element(const int& index)
{
    // If index is outside boundary then throw an error
    if (index < 0 || index >= elements()) {
        throw GException::out_of_range(G_ELEMENT1, index, 0, elements()-1);
    }

    // Get the requested child element
    GXmlElement* element  = NULL;
    int          elements = 0;
    for (int i = 0; i < m_nodes.size(); ++i) {
        GXmlElement* src = dynamic_cast<GXmlElement*>(m_nodes[i]);
        if (src != NULL) {
            if (elements == index) {
                element = src;
                break;
            }
            elements++;
        }
    }

    // Return child element
    return element;
}


/***********************************************************************//**
 * @brief Return pointer to GXMLElement child (const variant)
 *
 * @param[in] index Node index [0,...,elements()-1].
 * @return Pointer to GXMLElement child.
 *
 * @exception GException::out_of_range
 *            Child element index is out of range.
 *
 * Returns a pointer to the child number @p index of the XML node. An
 * exception will be thrown if the @p index is not valid.
 ***************************************************************************/
const GXmlElement* GXmlNode::element(const int& index) const
{
    // If index is outside boundary then throw an error
    if (index < 0 || index >= elements()) {
        throw GException::out_of_range(G_ELEMENT1, index, 0, elements()-1);
    }

    // Get the requested child element
    const GXmlElement* element  = NULL;
    int                elements = 0;
    for (int i = 0; i < m_nodes.size(); ++i) {
        const GXmlElement* src = dynamic_cast<const GXmlElement*>(m_nodes[i]);
        if (src != NULL) {
            if (elements == index) {
                element = src;
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
 * @param[in] index Node index [0,...,elements()-1].
 *
 * @exception GException::xml_name_not_found
 *            Child element name not found.
 * @exception GException::out_of_range
 *            Child element index is out of range.
 *
 * Returns a pointer to the child number @p index with @p name of the XML
 * node. An exception will be thrown if the @p index is not valid.
 ***************************************************************************/
GXmlElement* GXmlNode::element(const std::string& name, const int& index)
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

    // Get the requested child element
    GXmlElement* element  = NULL;
    int          elements = 0;
    for (int i = 0; i < m_nodes.size(); ++i) {
        GXmlElement* src = dynamic_cast<GXmlElement*>(m_nodes[i]);
        if (src != NULL) {
            if (src->name() == name) {
                if (elements == index) {
                    element = src;
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
 * @brief Return pointer on GXMLElement child of a given name (const variant)
 *
 * @param[in] name Name of child element.
 * @param[in] index Node index [0,...,elements()-1].
 *
 * @exception GException::xml_name_not_found
 *            Child element name not found.
 * @exception GException::out_of_range
 *            Child element index is out of range.
 *
 * Returns a pointer to the child number @p index with @p name of the XML
 * node. An exception will be thrown if the @p index is not valid.
 ***************************************************************************/
const GXmlElement* GXmlNode::element(const std::string& name, const int& index) const
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

    // Get the requested child element
    const GXmlElement* element  = NULL;
    int                elements = 0;
    for (int i = 0; i < m_nodes.size(); ++i) {
        const GXmlElement* src = dynamic_cast<const GXmlElement*>(m_nodes[i]);
        if (src != NULL) {
            if (src->name() == name) {
                if (elements == index) {
                    element = src;
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
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing XML information.
 ***************************************************************************/
std::string GXmlNode::print(const GChatter& chatter) const
{
    // Set result string
    std::string result = print(chatter, 0);

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
