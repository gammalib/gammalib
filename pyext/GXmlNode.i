/***************************************************************************
 *                 GXmlNode.i - Abstract XML node base class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GXmlNode.i
 * @brief Abstract XML node base class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXmlNode.hpp"
#include "GXmlElement.hpp"
#include "GXmlComment.hpp"
#include "GXmlDocument.hpp"
#include "GXmlPI.hpp"
#include "GXmlText.hpp"
%}

/* __ Typemaps ___________________________________________________________ */
%typemap(out) GXmlNode* {
    if (dynamic_cast<GXmlElement*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GXmlElement, 0 |  0 );
    }
    else if (dynamic_cast<GXmlComment*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GXmlComment, 0 |  0 );
    }
    else if (dynamic_cast<GXmlDocument*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GXmlDocument, 0 |  0 );
    }
    else if (dynamic_cast<GXmlPI*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GXmlPI, 0 |  0 );
    }
    else if (dynamic_cast<GXmlText*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GXmlText, 0 |  0 );
    }
    else {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GXmlNode, 0 |  0 );
    }
}


/***********************************************************************//**
 * @class GXmlNode
 *
 * @brief Abstract XML node base class
 ***************************************************************************/
class GXmlNode : public GBase {
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
    virtual void      clear(void) = 0;
    virtual GXmlNode* clone(void) const = 0;
    virtual void      write(GUrl& url, const int& indent) const = 0;
    virtual NodeType  type(void) const = 0;
    
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
    char *__str__() {
        return tochar(self->print());
    }
};
