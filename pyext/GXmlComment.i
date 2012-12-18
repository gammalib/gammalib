/***************************************************************************
 *                  GXmlComment.i - XML comment node class                 *
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
 * @file GXmlComment.i
 * @brief XML comment class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXmlComment.hpp"
#include "GTools.hpp"
%}

/* __ Includes ___________________________________________________________ */
%include "GTypemaps.i"


/***********************************************************************//**
 * @class GXmlComment
 *
 * @brief XML comment node class
 *
 * This class implements a XML comment. The comment text is stored without
 * the <!-- --> brackets.
 ***************************************************************************/
class GXmlComment : public GXmlNode {

public:
    // Constructors and destructors
    GXmlComment(void);
    GXmlComment(const GXmlComment& node);
    explicit GXmlComment(const std::string& segment);
    virtual ~GXmlComment(void);

    // Implemented virtual methods
    virtual void         clear(void);
    virtual GXmlComment* clone(void) const;
    virtual void         write(FILE* fptr, int indent = 0) const;
    virtual NodeType     type(void) const;

    // Other methods
    const std::string& comment(void) const;
    void               comment(const std::string& comment);
};


/***********************************************************************//**
 * @brief GXmlComment class extension
 ***************************************************************************/
%extend GXmlComment {
    char *__str__() {
        return tochar(self->print());
    }
};
