/***************************************************************************
 *                          GXml.hpp - XML class                           *
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
 * @file GXml.hpp
 * @brief XML class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GXML_HPP
#define GXML_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GContainer.hpp"
#include "GUrl.hpp"
#include "GXmlNode.hpp"
#include "GXmlDocument.hpp"
#include "GXmlElement.hpp"
#include "GXmlText.hpp"

/* __ Forward declarations _______________________________________________ */


/***********************************************************************//**
 * @class GXml
 *
 * @brief XML class
 *
 * This class holds the content of an Extensible Markup Language (XML)
 * document. An XML document is composed of a list of nodes, each of which
 * may contain lists of nodes, which again can contain list of nodes and so
 * on. This produces a tree made of nodes with an arbitrary complexity.
 * Nodes that do not contain any further nodes are the endpoints of the tree
 * that are called leafs.
 *
 * An example of an XML document is shown below
 *
 *     <?xml version="1.0" encoding="UTF-8" ?>
 *     <element type="Measurement">
 *       <parameter name="Flux" value="1.0"/>
 *     </element>
 *     <element>
 *       <list>
 *         <string>This is a text</string>
 *         <integer>17</integer>
 *       </list>
 *     </element>
 *
 * An XML document is a plain ASCII file. Every XML document begins by a
 * declaration of the XML version and an optional information about the
 * encoding of the text.
 *
 * The XML document is structured using @b tags. A tag is a markup construct
 * that begins with @< and ends with @>. Tags come in three flavors:
 * - start-tags; for example: @<section@>
 * - end-tags; for example: @</section@>
 * - empty-element tags; for example: @<line-break /@>
 * 
 * The header line of an XML document is an empty-element tag. Each tag may
 * contain an arbitrary number of @b attributes of the form
 *
 *     version="1.0"
 *
 * Alternative quotes ' are also allowed.
 *
 * A logical document component which either begins with a start-tag and ends
 * with a matching end-tag or consists only of an empty-element tag is called
 * an @b element. The characters between the start- and end-tags, if any, are
 * the element's content, and may contain markup, including other elements,
 * which are called @b child @b elements. An example of an element is
 *
 *  
 *     <string>This is a text</string>
 *
 * Another is
 *
 *     <parameter name="Flux" value="1.0"/>
 *
 * This last example has two attributes. The first word in an element is
 * called the @b element @b name. Every element has a name. Elements can
 * therefore be accessed by name, however, several elements with the same
 * name may exist.
 *
 * GammaLib implements the XML document in form of a master class GXml
 * which contains the root node of the document. The root node is
 * realized by the GXmlDocument class. GXmlDocument derives from GXmlNode,
 * which defines the abstract interface for all XML nodes. The following
 * XML nodes exist:
 * - GXmlElement: implements an element
 * - GXmlText: implements a text leaf
 * - GXmlPI: implements a Processing Instruction
 * - GXmlComment: implement a comment
 *
 * XML element attributes are implemented using GXmlAttribute.
 *
 * The GXml class provides methods to access the nodes of the XML document,
 * and to load and save the document from a URL. GXml derives from GContainer
 * as it behaves like a container class. It does, however, not contain an
 * explicit list, as the only data member of the class is a single instance
 * of GXmlDocument. GXmlDocument contains the hierachical list of all XML
 * nodes.
 *
 * To manipulate the child nodes of GXmlDocument, the usual container class
 * methods are available. size() returns the number of child nodes that
 * exist in the XML document root (in the above example there would be two
 * child elements with name @p element). The child nodes are accessed using
 * the operator[]. The is_empty() method checks whether the document root has
 * no children.
 *
 * The set() method allows to set a specific child node, the append() method
 * appends a child node to the XML document root. There is a second variant
 * of the append() method that appends a element of type GXmlElement to the
 * XML document root and that returns a pointer to this element. As argument,
 * this second method takes a text string that defines the element name and
 * attributes. For example, the first node in the above example could have
 * been generated using
 *
 *     GXml xml;
 *     xml.append("element type=\"Measurement\"");
 *
 * Note that the @< and @> symbols are not part of the text string that is
 * passed to the append() method.
 *
 * The insert() method inserts a child node at the specified index in the
 * XML root document. The remove() method removes the child node at the
 * specified index from the document root. The reserve() method reserves
 * space for a specified number of child nodes in the XML document root.
 * And the extend() method appends all child nodes that are found in the
 * specified argument to the XML document root.
 *
 * Most of the nodes encountered in a XML document will be @b element
 * @b nodes, hence special methods for handling element nodes have been
 * implemented. The elements() method returns the number of child elements
 * that are present in the XML document root. A variant that takes a string
 * argument counts the number of child elements with a given name. The
 * element() method allows to access the child elements. There are again
 * two flavours, one that simply takes an index to loop over all child
 * elements, and another that takes a name and a index to loop over all
 * child elements of a given name. For element access, non-const and const
 * variants exist.
 *
 * Finally, the load() and save() methods enable loading and saving a XML
 * document from and to disk. The read() and write() method, which do the
 * actual job of reading and writing, operate on general Unified Resource
 * Locators, so that XML documents can also through other media than flat
 * files.
 *
 * The print() method does not print the full XML document into a string
 * but shows a concise summary that reflects the tree structure. For writing
 * the document in a string, use the write() method to write in a string URL
 * object (GUrlString).
 ***************************************************************************/
class GXml : public GContainer {

    // Friend classes
    friend class GXmlNode;
    friend class GXmlDocument;
    friend class GXmlText;

public:
    // Constructors and destructors
    GXml(void);
    GXml(const GXml& xml);
    explicit GXml(const std::string& xml);
    virtual ~GXml(void);

    // Operators
    GXml&           operator=(const GXml& xml);
    GXmlNode*       operator[](const int& index);
    const GXmlNode* operator[](const int& index) const;

    // Methods
    void               clear(void);
    GXml*              clone(void) const;
    int                size(void) const { return m_root.size(); }
    bool               is_empty(void) const { return m_root.is_empty(); }
    GXmlNode*          set(const int& index, const GXmlNode& node);
    GXmlNode*          append(const GXmlNode& node);
    GXmlElement*       append(const std::string& segment);
    GXmlNode*          insert(const int& index, const GXmlNode& node);
    void               remove(const int& index);
    void               reserve(const int& num);
    void               extend(const GXmlNode& node);
    int                elements(void) const;
    int                elements(const std::string& name) const;
    GXmlElement*       element(const int& index);
    const GXmlElement* element(const int& index) const;
    GXmlElement*       element(const std::string& name, const int& index);
    const GXmlElement* element(const std::string& name, const int& index) const;
    void               load(const std::string& filename);
    void               save(const std::string& filename);
    void               read(const GUrl& url);
    void               write(GUrl& url, const int& indent = 0) const;
    std::string        print(const GChatter& chatter = NORMAL) const;
    std::string        print(const GChatter& chatter = NORMAL,
                             const int&      indent = 0) const;

protected:
    // Protected enumerators
    enum MarkupType {
        MT_ELEMENT_START,
        MT_ELEMENT_END,
        MT_ELEMENT_EMPTY,
        MT_COMMENT,
        MT_DECLARATION,
        MT_PROCESSING,
        MT_INVALID
    };

    // Protected methods
    void       init_members(void);
    void       copy_members(const GXml& xml);
    void       free_members(void);
    void       parse(const GUrl& url);
    void       process_markup(GXmlNode** current, const std::string& segment);
    void       process_text(GXmlNode** current, const std::string& segment);
    MarkupType get_markuptype(const std::string& segment) const;

    // Protected members
    GXmlDocument m_root;   //!< Root document node
};

#endif /* GXML_HPP */
