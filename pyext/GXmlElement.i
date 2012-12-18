/***************************************************************************
 *             GXmlElement.i - XML element node class definition           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
 * @file GXmlElement.i
 * @brief XML element node class Python interface definition
  * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXmlElement.hpp"
#include "GTools.hpp"
%}

/* __ Includes ___________________________________________________________ */
%include "GTypemaps.i"


/***********************************************************************//**
 * @class GXmlElement
 *
 * @brief XML element node class
 *
 * This class implements an XML element with it's associated attributes.
 ***************************************************************************/
class GXmlElement : public GXmlNode {

public:
    // Constructors and destructors
    GXmlElement(void);
    GXmlElement(const GXmlElement& node);
    explicit GXmlElement(const std::string& segment);
    virtual ~GXmlElement(void);

    // Implemented virtual methods
    virtual void         clear(void);
    virtual GXmlElement* clone(void) const;
    virtual void         write(FILE* fptr, int indent = 0) const;
    virtual NodeType     type(void) const;

    // Methods
    std::string name(void) const;
    std::string attribute(const std::string& name) const;
    GXmlNode*   parent(void) const;
    void        name(const std::string& name);
    void        parent(GXmlNode* node);
    void        attribute(const std::string& name, const std::string& value);
};


/***********************************************************************//**
 * @brief GXmlElement class extension
 ***************************************************************************/
%extend GXmlElement {
    char *__str__() {
        return tochar(self->print());
    }
};
