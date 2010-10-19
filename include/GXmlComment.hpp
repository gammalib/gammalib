/***************************************************************************
 *           GXmlComment.hpp - XML comment node class definition           *
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
 * @file GXmlComment.hpp
 * @brief XML comment node class definition
 * @author J. Knodlseder
 */

#ifndef GXMLCOMMENT_HPP
#define GXMLCOMMENT_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GXmlNode.hpp"


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

    // Operators
    GXmlComment& operator= (const GXmlComment& node);

    // Implemented virtual methods
    void     clear(void);
    void     write(FILE* fptr, int indent = 0) const;
    void     print(std::ostream& os, int indent = 0) const;
    NodeType type(void) const { return NT_COMMENT; }

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GXmlComment& node);
    void         free_members(void);
    GXmlComment* clone(void) const;
    void         parse(const std::string& segment);
    std::string  comment(void) const { return m_comment; }
    void         comment(const std::string& comment) { m_comment=comment; }

    // Protected data members
    std::string m_comment;       //!< Comment (excluding brackets)
};

#endif /* GXMLCOMMENT_HPP */
