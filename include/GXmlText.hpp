/***************************************************************************
 *              GXmlText.hpp - XML text node class definition              *
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
 * @file GXmlText.hpp
 * @brief XML text node class definition
 * @author J. Knodlseder
 */

#ifndef GXMLTEXT_HPP
#define GXMLTEXT_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GXmlNode.hpp"


/***********************************************************************//**
 * @class GXmlText
 *
 * @brief XML text node class interface defintion.
 *
 * This class implements a text segment of a XML document.
 ***************************************************************************/
class GXmlText : public GXmlNode {

public:
    // Constructors and destructors
    GXmlText(void);
    GXmlText(const GXmlText& node);
    GXmlText(const std::string& text);
    ~GXmlText(void);

    // Operators
    GXmlText& operator= (const GXmlText& node);

    // Implemented virtual methods
    void     clear(void);
    void     write(FILE* fptr, int indent = 0) const;
    void     print(std::ostream& os, int indent = 0) const;
    NodeType type(void) const { return NT_TEXT; }

protected:
    // Protected methods
    void      init_members(void);
    void      copy_members(const GXmlText& node);
    void      free_members(void);
    GXmlText* GXmlText::clone(void) const;

    // Protected data members
    std::string m_text;       //!< Text
};

#endif /* GXMLTEXT_HPP */
