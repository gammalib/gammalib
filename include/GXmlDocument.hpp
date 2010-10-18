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
#include <string>
#include "GXmlNode.hpp"
#include "GXmlAttribute.hpp"


/***********************************************************************//**
 * @class GXmlDocument
 *
 * @brief XML document node class interface defintion.
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

    // Operators
    GXmlDocument& operator= (const GXmlDocument& node);

    // Implemented virtual methods
    void     clear(void);
    void     write(FILE* fptr, int indent = 0) const;
    void     print(std::ostream& os, int indent = 0) const;
    NodeType type(void) const { return NT_DOCUMENT; }

    // Methods
    std::string version(void) const { return m_version.value(); }
    std::string encoding(void) const { return m_encoding.value(); }
    std::string standalone(void) const { return m_standalone.value(); }
    void        version(const std::string& version) { m_version.value(version); }
    void        encoding(const std::string& encoding) { m_encoding.value(encoding); }
    void        standalone(const std::string& standalone) { m_standalone.value(standalone); }

protected:
    // Protected methods
    void          init_members(void);
    void          copy_members(const GXmlDocument& node);
    void          free_members(void);
    GXmlDocument* GXmlDocument::clone(void) const;

    // Protected data members
    GXmlAttribute m_version;      //!< XML version ("1.0", "1.1")
    GXmlAttribute m_encoding;     //!< Encoding (e.g. "UTF-8")
    GXmlAttribute m_standalone;   //!< Standalone ("yes", "no") 
};

#endif /* GXMLDOCUMENT_HPP */
