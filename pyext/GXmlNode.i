/***************************************************************************
 *                 GXmlNode.i - XML node base class definition             *
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
 * @file GXmlNode.i
 * @brief GXmlNode class python bindings
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXmlNode.hpp"
%}


/***********************************************************************//**
 * @class GXmlNode
 *
 * @brief XML node base class interface defintion.
 *
 * This class defines an abstract node of a XML document.
 ***************************************************************************/
class GXmlNode {
public:
    // Constructors and destructors
    GXmlNode(void);
    GXmlNode(const GXmlNode& node);
    virtual ~GXmlNode(void);

    // Public enumerators
    enum NodeType {
        NT_DOCUMENT,
        NT_ELEMENT,
        NT_COMMENT,
        NT_UNKNOWN,
        NT_TEXT,
        NT_DECLARATION,
        NT_PI,
        NT_TYPECOUNT
    };

    // Pure virtual methods
    virtual void     clear(void) = 0;
    virtual void     write(FILE* fptr, int indent = 0) const = 0;
    //virtual void     print(std::ostream& os, int indent = 0) const = 0;
    virtual NodeType type(void) const = 0;
    
    // Methods
    void      append(GXmlNode* node);
    int       children(void) const;
    GXmlNode* child(int index) const;
    int       elements(void) const;
    int       elements(const std::string& name) const;
    GXmlNode* element(int index) const;
    GXmlNode* element(const std::string& name, int index) const;
};


/***********************************************************************//**
 * @brief GXmlNode class extension
 ***************************************************************************/
%extend GXmlNode {
//    char *__str__() {
//        static std::string result = self->print();
//        return ((char*)result.c_str());
//    }
};
