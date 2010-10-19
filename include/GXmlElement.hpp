/***************************************************************************
 *            GXmlElement.hpp - XML element node class definition          *
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
 * @file GXmlElement.hpp
 * @brief XML element node class definition
 * @author J. Knodlseder
 */

#ifndef GXMLELEMENT_HPP
#define GXMLELEMENT_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include <iostream>
#include "GXml.hpp"
#include "GXmlNode.hpp"
#include "GXmlAttribute.hpp"


/***********************************************************************//**
 * @class GXmlElement
 *
 * @brief XML element node class interface defintion.
 *
 * This class implements an XML element with it's associated attributes.
 ***************************************************************************/
class GXmlElement : public GXmlNode {

    // Friend classes
    friend class GXml;

public:
    // Constructors and destructors
    GXmlElement(void);
    GXmlElement(const GXmlElement& node);
    GXmlElement(const std::string& segment);
    ~GXmlElement(void);

    // Operators
    GXmlElement& operator= (const GXmlElement& node);

    // Implemented virtual methods
    void     clear(void);
    void     write(FILE* fptr, int indent = 0) const;
    void     print(std::ostream& os, int indent = 0) const;
    NodeType type(void) const { return NT_ELEMENT; }

    // Methods
    std::string name(void) const { return m_name; }
    std::string attribute(const std::string& name) const;
    GXmlNode*   parent(void) const { return m_parent; }
    void        name(const std::string& name) { m_name=name; }
    void        parent(GXmlNode* node) { m_parent = node; }
    void        attribute(const std::string& name, const std::string& value);

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GXmlElement& node);
    void         free_members(void);
    GXmlElement* clone(void) const;
    void         parse_start(const std::string& segment);
    void         parse_stop(const std::string& segment);
    void         parse_attribute(size_t* pos, const std::string& segment);

    // Protected data members
    std::string                 m_name;       //!< Element name
    GXmlNode*                   m_parent;     //!< Pointer on parent node
    std::vector<GXmlAttribute*> m_attr;       //!< Attributes
};

#endif /* GXMLELEMENT_HPP */
