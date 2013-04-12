/***************************************************************************
 *          GXmlElement.cpp - XML element node class implementation        *
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
 * @file GXmlElement.cpp
 * @brief XML element node class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GXmlElement.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_PARSE_START                "GXmlElement::parse_start(std::string&)"
#define G_PARSE_STOP                  "GXmlElement::parse_stop(std::string&)"
#define G_PARSE_ATTRIBUTE            "GXmlElement::parse_attribute(size_t*, "\
                                                              "std::string&)"

/* __ Constants __________________________________________________________ */
const int g_indent = 2;                      //!< Indent for XML file writing

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
GXmlElement::GXmlElement(void) : GXmlNode()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] node XML element.
 ***************************************************************************/
GXmlElement::GXmlElement(const GXmlElement& node) : GXmlNode(node)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(node);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Segment constructor
 *
 * @param[in] segment XML segment.
 *
 * Constructs a XML element from a text @p segment. The text segment is
 * parsed and the element name and attributes are extracted using the
 * parse_start() method.
 ***************************************************************************/
GXmlElement::GXmlElement(const std::string& segment) : GXmlNode()
{
    // Initialise members
    init_members();

    // Parse start element
    parse_start(segment);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GXmlElement::~GXmlElement(void)
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
 * @param[in] node XML element.
 * @return XML element.
 ***************************************************************************/
GXmlElement& GXmlElement::operator=(const GXmlElement& node)
{
    // Execute only if object is not identical
    if (this != &node) {

        // Copy base class members
        this->GXmlNode::operator=(node);

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
 * @brief Clear XML element
 *
 * Resets the XML element to a clean initial state.
 ***************************************************************************/
void GXmlElement::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GXmlNode::free_members();

    // Initialise members
    this->GXmlNode::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone XML element
 *
 * @return Pointer to deep copy of XML element.
 ***************************************************************************/
GXmlElement* GXmlElement::clone(void) const
{
    // Clone element
    return new GXmlElement(*this);
}


/***********************************************************************//**
 * @brief Return attribute value
 *
 * @param[in] name Attribute name.
 * @return String containing attribute value.
 *
 * Returns the value of the attribute @p name. If the requested attribute was
 * not found an empty string is returned.
 ***************************************************************************/
std::string GXmlElement::attribute(const std::string& name) const
{
    // Initialise empty value (i.e. attribute not found)
    std::string value = "";

    // Search attribute value in list of attributes
    for (int i = 0; i < m_attr.size(); ++i) {
        if (m_attr[i]->name() == name) {
            value = m_attr[i]->value();
            break;
        }
    }

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Set attribute value
 *
 * @param[in] name Attribute name.
 * @param[in] value Attribute value.
 *
 * Sets an attribute of the element. If the attribute name exists the value
 * is modified. If the attribute does not yet exist it is created and
 * added to the list of attributes.
 *
 * Note that this logical assures that only one attribute with a given name
 * will exist in the element.
 ***************************************************************************/
void GXmlElement::attribute(const std::string& name, const std::string& value)
{
    // Initialise attribute NULL pointer
    GXmlAttribute* attr = NULL;

    // Search attribute name in list of attributes
    for (int i = 0; i < m_attr.size(); ++i) {
        if (m_attr[i]->name() == name) {
            attr = m_attr[i];
            break;
        }
    }

    // If no attribute with specified name has been found then add a new
    // attribute to the list of attributes
    if (attr == NULL) {
        attr = new GXmlAttribute;
        attr->name(name);
        m_attr.push_back(attr);
    }

    // Set or update value of attribute
    attr->value(value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Remove attribute from element
 *
 * @param[in] name Attribute name.
 *
 * Remove the attribute with @p name from the XML element. If the requested
 * attribute was not found the method does nothing.
 ***************************************************************************/
void GXmlElement::remove_attribute(const std::string& name)
{
    // Do nothing if there are no attributes
    if (!m_attr.empty()) {

        // Store number of attributes.
        int num = m_attr.size();

        // Search attribute name in list of attributes and erase attribute
        // when it has been found. Note that if several attributes with the
        // same name exist (which should never be the case!), only the
        // first attribute is removed
        for (int i = 0; i < num; ++i) {
            if (m_attr[i]->name() == name) {
                m_attr.erase(m_attr.begin() + i);
                break;
            }
        }

    } // endif: there were attributes

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write element into URL
 *
 * @param[in] url Unified Resource Locator.
 * @param[in] indent Text indentation (default = 0).
 *
 * Writes the element into a Unified Resource Locator.
 ***************************************************************************/
void GXmlElement::write(GUrl& url, const int& indent) const
{
    // Prepend indentation
    for (int k = 0; k < indent; ++k) {
        url.printf(" ");
    }

    // Write element name into URL
    url.printf("<%s", m_name.c_str());

    // Write attributes into URL
    for (int k = 0; k < m_attr.size(); ++k) {
        m_attr[k]->write(url);
    }

    // If there are no children then write an empty tag
    if (isempty()) {
        url.printf(" />\n");
    }

    // ... otherwise finish start tag, write children and write end tag
    else {

        // Case A: The element contains a single text leaf
        if ((m_nodes.size() == 1) && (m_nodes[0]->type() == NT_TEXT)) {

            // Finish start tag
            url.printf(">");

            // Write text leaf
            m_nodes[0]->write(url, 0);

            // Write end tag
            url.printf("</%s>\n", m_name.c_str());
        }
        
        // Case B: ... otherwise it contains markup
        else {

            // Finish start tag
            url.printf(">\n");

            // Write children in file
            for (int i = 0; i < m_nodes.size(); ++i) {
                m_nodes[i]->write(url, indent+g_indent);
                if (m_nodes[i]->type() == NT_TEXT) {
                    url.printf("\n");
                }
            }

            // Write end tag
            for (int k = 0; k < indent; ++k) {
                url.printf(" ");
            }
            url.printf("</%s>\n", m_name.c_str());
        
        } // endelse: element contained markup
        
    } // endelse: finished start tag

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print XML element
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @param[in] indent Text indentation (defaults to 0).
 * @return String containing XML element
 ***************************************************************************/
std::string GXmlElement::print(const GChatter& chatter,
                               const int&      indent) const
{
    // Allocate result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Initialise result string
        result = gammalib::fill(" ", indent);

        // Append element to string
        result.append("GXmlElement::"+m_name);
        for (int k = 0; k < m_attr.size(); ++k) {
            result.append(m_attr[k]->print(chatter));
        }

        // Append children
        for (int i = 0; i < m_nodes.size(); ++i) {
            result.append("\n" + m_nodes[i]->print(chatter, indent+g_indent));
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
void GXmlElement::init_members(void)
{
    // Initialise members
    m_name.clear();
    m_attr.clear();
    m_parent = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] node XML element.
 *
 * This method copies all class members. XML attributes are cloned.
 *
 * @todo Is copying the parent correct?
 ***************************************************************************/
void GXmlElement::copy_members(const GXmlElement& node)
{
    // Copy members
    m_name   = node.m_name;
    m_parent = node.m_parent;

    // Copy attribute container
    m_attr.clear();
    for (int i = 0; i < node.m_attr.size(); ++i) {
        m_attr.push_back((node.m_attr[i]->clone()));
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
void GXmlElement::free_members(void)
{
    // Free attributes
    for (int i = 0; i < m_attr.size(); ++i) {
        delete m_attr[i];
        m_attr[i] = NULL;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Parse element start segment string
 *
 * @param[in] segment Segment string.
 *
 * @exception GException::xml_syntax_error
 *            XML syntax error.
 *
 * Parse the segment string and set class members based on the information
 * that is found. The method also performs syntax checking. It does not
 * require brackets to be set.
 ***************************************************************************/
void GXmlElement::parse_start(const std::string& segment)
{
    // Initialize position check
    std::size_t pos_start = 0;

    // Get length of segment
    int n = segment.length();

    // Throw an error is segment is empty
    if (n < 1) {
        throw GException::xml_syntax_error(G_PARSE_START, segment,
                          "no element name specified");
    }

    // If string starts with brackets then check that the brackets are
    // valid comment brackets
    if (segment[0] == '<') {
        if (n < 2 || (segment.compare(0,1,"<") != 0) ||
                     (segment.compare(n-1,1,">") != 0)) {
            throw GException::xml_syntax_error(G_PARSE_START, segment,
                                               "invalid tag brackets");
        }
        pos_start = 1;
    } // endif: there were brackets

    // Extract element name
    std::size_t pos = segment.find_first_of("\x20\x09\x0d\x0a>", 1);
    if (pos == pos_start) {
        throw GException::xml_syntax_error(G_PARSE_START, segment,
                          "no whitespace allowed before element name");
    }
    if (pos == std::string::npos) {
        if (pos_start == 1) {
            throw GException::xml_syntax_error(G_PARSE_START, segment,
                              "element name not found");
        }
    }
    m_name = segment.substr(pos_start, pos-pos_start);

    // Extract attributes
    while (pos != std::string::npos) {
        parse_attribute(&pos, segment);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Parse element stop segment string
 *
 * @param[in] segment Segment string.
 *
 * @exception GException::xml_syntax_error
 *            XML syntax error.
 *
 * Parse the stop segment string and verify the syntax.
 ***************************************************************************/
void GXmlElement::parse_stop(const std::string& segment)
{
    // Get length of segment
    int n = segment.length();

    // Check on existence of brackets
    if (n < 3 || (segment.compare(0,2,"</") != 0) ||
                 (segment.compare(n-1,1,">") != 0)) {
        throw GException::xml_syntax_error(G_PARSE_STOP, segment,
                          "incorrect or missing tag brackets");
    }

    // Extract and verify element name
    size_t pos = segment.find_first_of("\x20\x09\x0d\x0a>", 2);
    if (pos == 2) {
        throw GException::xml_syntax_error(G_PARSE_STOP, segment,
                          "no whitespace allowed after \"</\"");
    }
    if (pos == std::string::npos) {
        throw GException::xml_syntax_error(G_PARSE_STOP, segment,
                          "element name not found");
    }
    std::string name = segment.substr(2, pos-2);
    if (name != m_name) {
        throw GException::xml_syntax_error(G_PARSE_STOP, segment,
                          "invalid name in element stop tag"
                          " (found \""+name+"\", expected \""+m_name+"\"");
    }

    // Verify that no further characters exist in element stop tag
    size_t pos2 = segment.find_first_of("\x20\x09\x0d\x0a>", pos);
    if (pos2 != n-1) {
        throw GException::xml_syntax_error(G_PARSE_STOP, segment,
                          "invalid characters found after element name");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Parse element attribute
 *
 * @param[in] pos Start position in string.
 * @param[in] segment Segment string.
 *
 * @exception GException::xml_syntax_error
 *            XML syntax error.
 *
 * Parse the segment string for one attribute, and if attribute was found,
 * attach it to element.
 *
 * @todo Verify XML validity of attribute name and value
 ***************************************************************************/
void GXmlElement::parse_attribute(size_t* pos, const std::string& segment)
{
    // Main loop
    do {
        // Get substring for error message
        std::string error = segment.substr(*pos, segment.length()-*pos);

        // Find first character of name substring
        std::size_t pos_name_start = segment.find_first_not_of("\x20\x09\x0d\x0a/>?", *pos);
        if (pos_name_start == std::string::npos) {
            *pos = std::string::npos;
            continue;
        }

        // Find end of name substring
        std::size_t pos_name_end = segment.find_first_of("\x20\x09\x0d\x0a=", pos_name_start);
        if (pos_name_end == std::string::npos) {
            throw GException::xml_syntax_error(G_PARSE_ATTRIBUTE, error,
                              "invalid or missing attribute name");
        }

        // Find '=' character
        std::size_t pos_equal = segment.find_first_of("=", pos_name_end);
        if (pos_equal == std::string::npos) {
            throw GException::xml_syntax_error(G_PARSE_ATTRIBUTE, error,
                              "\"=\" sign not found for attribute");
        }

        // Find start of value substring
        std::size_t pos_value_start = segment.find_first_of("\x22\x27", pos_equal);
        if (pos_value_start == std::string::npos) {
            throw GException::xml_syntax_error(G_PARSE_ATTRIBUTE, error,
                              "invalid or missing attribute value start hyphen");
        }

        // Save hyphen character and step forward one character
        std::string hyphen = segment.substr(pos_value_start, 1);
        pos_value_start++;
        if (pos_value_start >= segment.length()) {
            throw GException::xml_syntax_error(G_PARSE_ATTRIBUTE, error,
                              "invalid or missing attribute value");
        }

        // Find end of value substring
        std::size_t pos_value_end = segment.find_first_of(hyphen, pos_value_start);
        if (pos_value_end == std::string::npos) {
            throw GException::xml_syntax_error(G_PARSE_ATTRIBUTE, error,
                              "invalid or missing attribute value end hyphen");
        }

        // Get name substring
        std::size_t n_name = pos_name_end - pos_name_start;
        if (n_name < 1) {
            throw GException::xml_syntax_error(G_PARSE_ATTRIBUTE, error,
                              "invalid or missing attribute name");
        }
        std::string name = segment.substr(pos_name_start, n_name);

        //@todo Check XML validity of attribute name

        // Get value substring length
        std::size_t n_value = pos_value_end - pos_value_start;
        //if (n_value < 0) {
        //    throw GException::xml_syntax_error(G_PARSE_ATTRIBUTE, error,
        //                      "invalid or missing attribute value");
        //}
        std::string value = segment.substr(pos_value_start-1, n_value+2);

        //@todo Check XML validity of attribute value

        // Allocate, set and append new attribute to element
        GXmlAttribute* attr  = new GXmlAttribute(name, value);
        m_attr.push_back(attr);

        // Update segment pointer
        pos_value_end++;
        if (pos_value_end >= segment.length()) {
            pos_value_end = std::string::npos;
        }
        *pos = pos_value_end;

    } while (0);

    // Return
    return;
}
