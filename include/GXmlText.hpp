/***************************************************************************
 *                   GXmlText.hpp - XML text node class                    *
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
 * @file GXmlText.hpp
 * @brief XML text node class interface definition
 * @author J. Knodlseder
 */

#ifndef GXMLTEXT_HPP
#define GXMLTEXT_HPP

/* __ Includes ___________________________________________________________ */
#include <cstdio>             // FILE*, std::fprintf
#include <string>
#include "GXmlNode.hpp"


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

    // Operators
    GXmlText& operator= (const GXmlText& node);

    // Implemented virtual methods
    virtual void        clear(void);
    virtual GXmlText*   clone(void) const;
    virtual void        write(FILE* fptr, int indent = 0) const;
    virtual std::string print(int indent = 0) const;
    virtual NodeType    type(void) const { return NT_TEXT; }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GXmlText& node);
    void free_members(void);

    // Protected data members
    std::string m_text;       //!< Text
};

#endif /* GXMLTEXT_HPP */
