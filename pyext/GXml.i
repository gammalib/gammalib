/***************************************************************************
 *                        GXml.i - XML class definition                    *
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
 * @file GXml.i
 * @brief GXml class python bindings
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXml.hpp"
%}


/***********************************************************************//**
 * @class GXml
 *
 * @brief XML class interface defintion.
 *
 * This class implements an XML document. It holds the root node.
 ***************************************************************************/
class GXml {
public:
    // Constructors and destructors
    GXml(void);
    GXml(const GXml& xml);
    GXml(const std::string& filename);
    ~GXml(void);

    // Methods
    void         clear(void);
    void         append(GXmlNode* node);
    void         load(const std::string& filename);
    void         save(const std::string& filename);
    int          children(void) const;
    GXmlNode*    child(int index) const;
    int          elements(void) const;
    int          elements(const std::string& name) const;
    GXmlElement* element(int index) const;
    GXmlElement* element(const std::string& name, int index) const;
};


/***********************************************************************//**
 * @brief GXml class extension
 ***************************************************************************/
%extend GXml {
//    char *__str__() {
//        static std::string result = self->print();
//        return ((char*)result.c_str());
//    }
};
