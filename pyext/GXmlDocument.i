/***************************************************************************
 *                 GXmlDocument.i - XML document node class                *
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
 * @file GXmlDocument.i
 * @brief XML document class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXmlDocument.hpp"
#include "GTools.hpp"
%}

/* __ Includes ___________________________________________________________ */
%include "GTypemaps.i"


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
    NodeType      type(void) const;

    // Methods
    std::string version(void) const;
    std::string encoding(void) const;
    std::string standalone(void) const;
    void        version(const std::string& version);
    void        encoding(const std::string& encoding);
    void        standalone(const std::string& standalone);
};


/***********************************************************************//**
 * @brief GXmlDocument class extension
 ***************************************************************************/
%extend GXmlDocument {
    char *__str__() {
        return tochar(self->print());
    }
};
