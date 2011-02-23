/***************************************************************************
 *                 GXmlDocument.i - XML document node class                *
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
 * @file GXmlDocument.i
 * @brief XML document class Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXmlDocument.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GXmlDocument
 *
 * @brief XML document node class
 *
 * This class implements the root node of an XML document. It contains the
 * three special attributes 'version', 'encoding', and 'standalone'.
 ***************************************************************************/
class GXmlDocument : public GXmlNode {

public:
    // Constructors and destructors
    GXmlDocument(void);
    GXmlDocument(const GXmlDocument& node);
    ~GXmlDocument(void);

    // Implemented virtual methods
    void          clear(void);
    GXmlDocument* clone(void) const;
    void          write(FILE* fptr, int indent = 0) const;
    //void          print(std::ostream& os, int indent = 0) const;
    NodeType      type(void) const;

    // Methods
    const std::string& version(void) const;
    const std::string& encoding(void) const;
    const std::string& standalone(void) const;
    void               version(const std::string& version);
    void               encoding(const std::string& encoding);
    void               standalone(const std::string& standalone);
};


/***********************************************************************//**
 * @brief GXmlDocument class extension
 ***************************************************************************/
%extend GXmlDocument {
//    char *__str__() {
//        return tochar(self->print());
//    }
};


/***********************************************************************//**
 * @brief GXmlDocument type casts
 ***************************************************************************/
%inline %{
    GXmlDocument* cast_GXmlDocument(GXmlNode* node) {
        if (node->type() != GXmlNode::NT_DOCUMENT)
            throw GException::xml_bad_node_type("cast_GXmlDocument(GXmlNode*)",
                                                "",
                                                "Expecting GXmlDocument node.");
        return dynamic_cast<GXmlDocument*>(node);
    }
%};
