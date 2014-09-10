/***************************************************************************
 *           GXmlAttribute.hpp - XML attribute class definition            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
 * @file GXmlAttribute.hpp
 * @brief XML attribute class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GXMLATTRIBUTE_HPP
#define GXMLATTRIBUTE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GUrl.hpp"


/***********************************************************************//**
 * @class GXmlAttribute
 *
 * @brief XML attribute class
 *
 * This class implements an attribute of an XML element. An attribute
 * consists of a name-value pair. There can only be one attribute with a
 * given name. Here an example of an XML element with attributes:
 *
 *     <element type="singleton" name="gamma-ray">
 *
 * The element has two attributes:
 * - @p type has a value of "singleton"
 * - @p name has a value of "gamma-ray"
 *
 * The hyphens are stored together with the attribute value. Allowed hyphens
 * are " and '.
 ***************************************************************************/
class GXmlAttribute : public GBase {

public:
    // Constructors and destructors
    GXmlAttribute(void);
    explicit GXmlAttribute(const GXmlAttribute& attr);
    explicit GXmlAttribute(const std::string& name, const std::string& value);
    virtual ~GXmlAttribute(void);

    // Operators
    GXmlAttribute& operator=(const GXmlAttribute& attr);

    // Methods
    void               clear(void);
    GXmlAttribute*     clone(void) const;
    std::string        classname(void) const;
    void               write(GUrl& url) const;
    const std::string& name(void) const { return m_name; }
    std::string        value(void) const;
    void               name(const std::string& name) { m_name=name; }
    void               value(std::string value);
    std::string        print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GXmlAttribute& attr);
    void free_members(void);

    // Protected data members
    std::string m_name;     //!< Attribute name
    std::string m_value;    //!< Attribute value including hyphens
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GXmlAttribute").
 ***************************************************************************/
inline
std::string GXmlAttribute::classname(void) const
{
    return ("GXmlAttribute");
}

#endif /* GXMLATTRIBUTE_HPP */
