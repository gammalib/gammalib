/***************************************************************************
 *            GXmlComment.i - XML comment node class definition            *
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
 * @file GXmlComment.i
 * @brief GXmlComment class python bindings
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXmlComment.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GXmlComment
 *
 * @brief XML comment node class interface defintion.
 *
 * This class implements a XML comment. The comment text is stored without
 * the <!-- --> brackets.
 ***************************************************************************/
class GXmlComment : public GXmlNode {

public:
    // Constructors and destructors
    GXmlComment(void);
    GXmlComment(const GXmlComment& node);
    GXmlComment(const std::string& segment);
    ~GXmlComment(void);

    // Implemented virtual methods
    void     clear(void);
    void     write(FILE* fptr, int indent = 0) const;
    //void     print(std::ostream& os, int indent = 0) const;
    NodeType type(void) const { return NT_COMMENT; }
};


/***********************************************************************//**
 * @brief GXmlComment class extension
 ***************************************************************************/
%extend GXmlComment {
//    char *__str__() {
//        return tochar(self->print());
//    }
    GXmlComment(const GXmlNode& node) {
        if (node.type() != GXmlNode::NT_COMMENT)
            throw GException::xml_bad_node_type("GXmlComment(GXmlNode&)",
                                                "",
                                                "Expecting GXmlComment node.");
        return (GXmlComment*)&node;
    }
};
