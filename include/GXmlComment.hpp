/***************************************************************************
 *           GXmlComment.hpp - XML comment node class definition           *
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
 * @file GXmlComment.hpp
 * @brief XML comment node class interface definition
 * @author J. Knodlseder
 */

#ifndef GXMLCOMMENT_HPP
#define GXMLCOMMENT_HPP

/* __ Includes ___________________________________________________________ */
#include <cstdio>             // FILE*, std::fprintf
#include <string>
#include "GXmlNode.hpp"


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

    // Operators
    GXmlComment& operator= (const GXmlComment& node);

    // Implemented virtual methods
    virtual void         clear(void);
    virtual GXmlComment* clone(void) const;
    virtual void         write(FILE* fptr, int indent = 0) const;
    virtual std::string  print(int indent = 0) const;
    virtual NodeType     type(void) const { return NT_COMMENT; }

    // Other methods
    const std::string& comment(void) const { return m_comment; }
    void               comment(const std::string& comment) { m_comment=comment; }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GXmlComment& node);
    void free_members(void);
    void parse(const std::string& segment);

    // Protected data members
    std::string m_comment;       //!< Comment (excluding brackets)
};

#endif /* GXMLCOMMENT_HPP */
