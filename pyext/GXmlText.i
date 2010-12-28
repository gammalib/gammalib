/***************************************************************************
 *               GXmlText.i - XML text node class definition               *
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
 * @file GXmlText.i
 * @brief GXmlText class python bindings
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXmlText.hpp"
%}


/***********************************************************************//**
 * @class GXmlText
 *
 * @brief XML text node class interface defintion.
 *
 * This class implements a text segment of a XML document.
 ***************************************************************************/
class GXmlText : public GXmlNode {
public:
    // Constructors and destructors
    GXmlText(void);
    GXmlText(const GXmlText& node);
    GXmlText(const std::string& text);
    ~GXmlText(void);

    // Implemented virtual methods
    void     clear(void);
    void     write(FILE* fptr, int indent = 0) const;
    //void     print(std::ostream& os, int indent = 0) const;
    NodeType type(void) const { return NT_TEXT; }
};


/***********************************************************************//**
 * @brief GXmlText class extension
 ***************************************************************************/
%extend GXmlText {
//    char *__str__() {
//        static std::string result = self->print();
//        return ((char*)result.c_str());
//    }
    GXmlText(const GXmlNode& node) {
        if (node.type() != GXmlNode::NT_TEXT)
            throw GException::xml_bad_node_type("GXmlText(GXmlNode&)",
                                                "",
                                                "Expecting GXmlText node.");
        return (GXmlText*)&node;
    }
};
