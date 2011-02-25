/***************************************************************************
 *                    GXmlText.i - XML text node class                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @file GXmlText.i
 * @brief XML text node class Python interface definition
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
 * @brief XML text node class
 *
 * This class implements a text segment of a XML document.
 ***************************************************************************/
class GXmlText : public GXmlNode {

public:
    // Constructors and destructors
    GXmlText(void);
    GXmlText(const GXmlText& node);
    explicit GXmlText(const std::string& text);
    virtual ~GXmlText(void);

    // Implemented virtual methods
    virtual void        clear(void);
    virtual GXmlText*   clone(void) const;
    virtual void        write(FILE* fptr, int indent = 0) const;
    virtual NodeType    type(void) const;
};


/***********************************************************************//**
 * @brief GXmlText class extension
 ***************************************************************************/
%extend GXmlText {
    char *__str__() {
        return tochar(self->print());
    }
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
