/***************************************************************************
 *                GXmlNode.hpp - XML node base class definition            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GXmlNode.hpp
 * @brief XML node base class definition
 * @author J. Knodlseder
 */

#ifndef GXMLNODE_HPP
#define GXMLNODE_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <iostream>

/***********************************************************************//**
 * @class GXmlNode
 *
 * @brief XML node base class interface defintion.
 ***************************************************************************/
class GXmlNode {

public:
    // Constructors and destructors
    GXmlNode(void);
    GXmlNode(const GXmlNode& node);
    ~GXmlNode(void);

    // Operators
    GXmlNode&       operator= (const GXmlNode& node);
    GXmlNode&       operator() (int index);
    const GXmlNode& operator() (int index) const;

    // Public enumerators
    enum NodeType {
        NT_DOCUMENT,
        NT_ELEMENT,
        NT_COMMENT,
        NT_UNKNOWN,
        NT_TEXT,
        NT_DECLARATION,
        NT_TYPECOUNT
    };

    // Pure virtual methods
    virtual void     clear(void) = 0;
    virtual void     print(std::ostream& os, int indent = 0) const = 0;
    virtual NodeType type(void) const = 0;
    
    // Methods
    void append(GXmlNode* node) { m_nodes.push_back(node); }
    int  size(void) const { return m_nodes.size(); }

protected:
    // Protected methods
    void              init_members(void);
    void              copy_members(const GXmlNode& node);
    void              free_members(void);
    virtual GXmlNode* GXmlNode::clone(void) const = 0;

    // Protected data members
    std::vector<GXmlNode*> m_nodes;    //!< Pointer to nodes contained in node
};

#endif /* GXMLNODE_HPP */
