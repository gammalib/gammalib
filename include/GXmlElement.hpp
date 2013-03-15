/***************************************************************************
 *            GXmlElement.hpp - XML element node class definition          *
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
 * @file GXmlElement.hpp
 * @brief XML element node class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GXMLELEMENT_HPP
#define GXMLELEMENT_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GUrl.hpp"
#include "GXmlNode.hpp"
#include "GXmlAttribute.hpp"


/***********************************************************************//**
 * @class GXmlElement
 *
 * @brief XML element node class
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
    explicit GXmlElement(const std::string& segment);
    virtual ~GXmlElement(void);

    // Operators
    GXmlElement& operator=(const GXmlElement& node);

    // Implemented virtual methods
    virtual void         clear(void);
    virtual GXmlElement* clone(void) const;
    virtual void         write(GUrl& url, const int& indent = 0) const;
    virtual NodeType     type(void) const { return NT_ELEMENT; }
    virtual std::string  print(const int& indent = 0) const;

    // Methods
    std::string  name(void) const { return m_name; }
    std::string  attribute(const std::string& name) const;
    GXmlNode*    parent(void) const { return m_parent; }
    void         name(const std::string& name) { m_name=name; }
    void         parent(GXmlNode* node) { m_parent = node; }
    void         attribute(const std::string& name, const std::string& value);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GXmlElement& node);
    void free_members(void);
    void parse_start(const std::string& segment);
    void parse_stop(const std::string& segment);
    void parse_attribute(size_t* pos, const std::string& segment);

    // Protected data members
    std::string                 m_name;       //!< Element name
    GXmlNode*                   m_parent;     //!< Pointer on parent node
    std::vector<GXmlAttribute*> m_attr;       //!< Attributes
};

#endif /* GXMLELEMENT_HPP */
