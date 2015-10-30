/***************************************************************************
 *                GXmlNode.hpp - Abstract XML node base class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2015 by Juergen Knoedlseder                         *
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
 * @file GXmlNode.hpp
 * @brief Abstract XML node base class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GXMLNODE_HPP
#define GXMLNODE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GContainer.hpp"
#include "GUrl.hpp"

/* __ Forward declarations _______________________________________________ */
class GXmlElement;


/***********************************************************************//**
 * @class GXmlNode
 *
 * @brief Abstract XML node base class
 *
 * This class defines the interface for all XML nodes. Each XML node can be
 * the container for a number of child nodes. The GXmlNode class is thus
 * designed as a container class that holds a list of pointers of type
 * GXmlNode. This allows implementing an arbitrary complex tree structure.
 *
 * The only member of this class is a list of XML node pointers. GXmlNode
 * handles the proper allocation and deallocation of the node memory.
 *
 * Most methods are identical to those of the GXml class. Please refer to
 * the documentation of this class for a description of the methods.
 ***************************************************************************/
class GXmlNode : public GContainer {

public:
    // Constructors and destructors
    GXmlNode(void);
    GXmlNode(const GXmlNode& node);
    virtual ~GXmlNode(void);

    // Operators
    GXmlNode&       operator=(const GXmlNode& node);
    GXmlNode*       operator[](const int& index);
    const GXmlNode* operator[](const int& index) const;

    // Public enumerators
    enum NodeType {
        NT_DOCUMENT,
        NT_ELEMENT,
        NT_COMMENT,
        NT_UNKNOWN,
        NT_TEXT,
        NT_DECLARATION,
        NT_PI,
        NT_TYPECOUNT
    };

    // Methods
    virtual void               clear(void) = 0;
    virtual GXmlNode*          clone(void) const = 0;
    virtual std::string        classname(void) const = 0;
    virtual int                size(void) const;
    virtual bool               is_empty(void) const;
    virtual GXmlNode*          set(const int& index, const GXmlNode& node);
    virtual GXmlNode*          append(const GXmlNode& node);
    virtual GXmlElement*       append(const std::string& segment);
    virtual GXmlNode*          insert(const int& index, const GXmlNode& node);
    virtual void               remove(const int& index);
    virtual void               reserve(const int& num);
    virtual void               extend(const GXmlNode& node);
    virtual int                elements(void) const;
    virtual int                elements(const std::string& name) const;
    virtual GXmlElement*       element(const int& index);
    virtual const GXmlElement* element(const int& index) const;
    virtual GXmlElement*       element(const std::string& name);
    virtual const GXmlElement* element(const std::string& name) const;
    virtual GXmlElement*       element(const std::string& name, const int& index);
    virtual const GXmlElement* element(const std::string& name, const int& index) const;
    virtual void               write(GUrl& url, const int& indent) const = 0;
    virtual NodeType           type(void) const = 0;
    virtual std::string        print(const GChatter& chatter = NORMAL,
                                     const int&      indent = 0) const = 0;
    virtual std::string        print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GXmlNode& node);
    void free_members(void);
    int  extract_index(std::string& tag) const;

    // Protected data members
    std::vector<GXmlNode*> m_nodes;   //!< Pointer to child nodes
};


/***********************************************************************//**
 * @brief Return number of child nodes
 *
 * @return Number of child nodes in node.
 ***************************************************************************/
inline
int GXmlNode::size(void) const
{
    return (int)m_nodes.size();
}


/***********************************************************************//**
 * @brief Signals if node has no child nodes
 *
 * @return True if node has no child nodes.
 ***************************************************************************/
inline
bool GXmlNode::is_empty(void) const
{
    return m_nodes.empty();
}


/***********************************************************************//**
 * @brief Reserve space for child nodes
 *
 * @param[in] num Number of child nodes for which space should be reserved.
 ***************************************************************************/
inline
void GXmlNode::reserve(const int& num)
{
    m_nodes.reserve(num);
    return;
}

#endif /* GXMLNODE_HPP */
