/***************************************************************************
 *        GXmlAttribute.hpp - XML attribute node class definition          *
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
 * @file GXmlAttribute.hpp
 * @brief XML attribute node class definition
 * @author J. Knodlseder
 */

#ifndef GXMLATTRIBUTE_HPP
#define GXMLATTRIBUTE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>


/***********************************************************************//**
 * @class GXmlAttribute
 *
 * @brief XML attribute node class interface defintion.
 *
 * This class implements the attribute of an XML elements. An attribute
 * consists of a name-value pair. Note that the hyphens are stored for the
 * attribute value. Allowed hyphens are " and '.
 ***************************************************************************/
class GXmlAttribute {

public:
    // Constructors and destructors
    GXmlAttribute(void);
    GXmlAttribute(const GXmlAttribute& attr);
    GXmlAttribute(const std::string& name, const std::string& value);
    ~GXmlAttribute(void);

    // Operators
    GXmlAttribute& operator= (const GXmlAttribute& attr);

    // Methods
    void        clear(void);
    void        write(FILE* fptr) const;
    void        print(std::ostream& os) const;
    std::string name(void) const { return m_name; }
    std::string value(void) const { return m_value; }
    void        name(const std::string& name) { m_name=name; }
    void        value(const std::string& value) { m_value=value; }

protected:
    // Protected methods
    void           init_members(void);
    void           copy_members(const GXmlAttribute& attr);
    void           free_members(void);
    GXmlAttribute* GXmlAttribute::clone(void) const;

    // Protected data members
    std::string m_name;       //!< Attribute name
    std::string m_value;      //!< Attribute value
};

#endif /* GXMLATTRIBUTE_HPP */
