/***************************************************************************
 *            GXmlElement.hpp - XML element node class definition          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2018 by Juergen Knoedlseder                         *
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
#include "GXmlNode.hpp"
#include "GXmlAttribute.hpp"

/* __ Forward declarations _______________________________________________ */
class GUrl;


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

    // Methods
    virtual void         clear(void);
    virtual GXmlElement* clone(void) const;
    virtual std::string  classname(void) const;
    const std::string&   name(void) const;
    void                 name(const std::string& name);
    int                  attributes(void) const;
    const GXmlAttribute* attribute(const int& index) const;
    std::string          attribute(const std::string& name) const;
    void                 attribute(const std::string& name,
                                   const std::string& value);
    bool                 has_attribute(const std::string& name) const;
    void                 remove_attribute(const std::string& name);
    virtual void         write(GUrl& url, const int& indent = 0) const;
    virtual NodeType     type(void) const;
    virtual std::string  print(const GChatter& chatter = NORMAL,
                               const int&     indent = 0) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GXmlElement& node);
    void free_members(void);
    void parse_start(const std::string& segment);
    void parse_stop(const std::string& segment);
    void parse_attribute(size_t* pos, const std::string& segment);

    // Protected data members
    std::string                 m_name;     //!< Element name
    std::vector<GXmlAttribute*> m_attr;     //!< Attributes
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GXmlElement").
 ***************************************************************************/
inline
std::string GXmlElement::classname(void) const
{
    return ("GXmlElement");
}


/***********************************************************************//**
 * @brief Return XML element name
 *
 * @return XML element name.
 ***************************************************************************/
inline
const std::string& GXmlElement::name(void) const
{
    return (m_name);
}


/***********************************************************************//**
 * @brief Set XML element name
 *
 * @param[in] name XML element name.
 ***************************************************************************/
inline
void GXmlElement::name(const std::string& name)
{
    m_name = name;
    return;
}


/***********************************************************************//**
 * @brief Return XML node type
 *
 * @return XML node type (NT_ELEMENT).
 ***************************************************************************/
inline
GXmlNode::NodeType GXmlElement::type(void) const
{
    return (NT_ELEMENT);
}


/***********************************************************************//**
 * @brief Return number of attributes
 *
 * @return Number of attributes.
 ***************************************************************************/
inline
int GXmlElement::attributes(void) const
{
    return (int)(m_attr.size());
}

#endif /* GXMLELEMENT_HPP */
