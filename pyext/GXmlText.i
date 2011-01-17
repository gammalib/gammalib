/***************************************************************************
 *               GXmlText.i - XML text node class definition               *
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
 * @file GXmlText.i
 * @brief GXmlText class python bindings
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXmlText.hpp"
#include "GTools.hpp"
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
//        return tochar(self->print());
//    }
};


/***********************************************************************//**
 * @brief GXmlText type casts
 ***************************************************************************/
%inline %{
    GXmlText* cast_GXmlText(GXmlNode* node) {
        if (node->type() != GXmlNode::NT_TEXT)
            throw GException::xml_bad_node_type("cast_GXmlText(GXmlNode*)",
                                                "",
                                                "Expecting GXmlText node.");
        return dynamic_cast<GXmlText*>(node);
    }
%};
