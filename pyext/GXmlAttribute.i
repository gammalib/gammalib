/***************************************************************************
 *         GXmlAttribute.i - XML attribute node class definition           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GXmlAttribute.i
 * @brief GXmlAttribute class python bindings
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXmlAttribute.hpp"
#include "GTools.hpp"
%}


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

    // Methods
    void        clear(void);
    void        write(FILE* fptr) const;
    //void        print(std::ostream& os) const;
    std::string name(void) const { return m_name; }
    std::string value(void) const;
    void        name(const std::string& name) { m_name=name; }
    void        value(std::string value);
};


/***********************************************************************//**
 * @brief GXmlAttribute class extension
 ***************************************************************************/
%extend GXmlAttribute {
//    char *__str__() {
//        return tochar(self->print());
//    }
};
