/***************************************************************************
 *           GXmlDocument.hpp - XML document node class definition         *
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
 * @file GXmlDocument.hpp
 * @brief XML document node class definition
 * @author J. Knodlseder
 */

#ifndef GXMLDOCUMENT_HPP
#define GXMLDOCUMENT_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GXmlNode.hpp"


/***********************************************************************//**
 * @class GXmlDocument
 *
 * @brief XML document node class interface defintion.
 ***************************************************************************/
class GXmlDocument : public GXmlNode {

public:
    // Constructors and destructors
    GXmlDocument(void);
    GXmlDocument(const GXmlDocument& node);
    ~GXmlDocument(void);

    // Operators
    GXmlDocument& operator= (const GXmlDocument& node);

    // Methods
    void     clear(void);
    void     print(std::ostream& os, int indent = 0) const;
    NodeType type(void) const { return NT_DOCUMENT; }

protected:
    // Protected methods
    void          init_members(void);
    void          copy_members(const GXmlDocument& node);
    void          free_members(void);
    GXmlDocument* GXmlDocument::clone(void) const;

    // Protected data members
};

#endif /* GXMLDOCUMENT_HPP */
