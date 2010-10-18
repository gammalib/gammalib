/***************************************************************************
 *                       GXml.hpp - XML class definition                   *
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
 * @file GXml.hpp
 * @brief XML class definition
 * @author J. Knodlseder
 */

#ifndef GXML_HPP
#define GXML_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GXmlNode.hpp"
#include "GXmlDocument.hpp"
#include "GXmlText.hpp"


/***********************************************************************//**
 * @class GXml
 *
 * @brief XML class interface defintion.
 *
 * This class implements an XML document. It holds the root node.
 ***************************************************************************/
class GXml {

    // Friend classes
    friend class GXmlNode;
    friend class GXmlDocument;
    friend class GXmlText;

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GXml& xml);

public:
    // Constructors and destructors
    GXml(void);
    GXml(const GXml& xml);
    GXml(const std::string& filename);
    ~GXml(void);

    // Operators
    GXml& operator= (const GXml& xml);

    // Methods
    void clear(void);
    void load(const std::string& filename);
    void save(const std::string& filename);

protected:
    // Protected enumerators
    enum MarkupType {
        MT_ELEMENT_START,
        MT_ELEMENT_END,
        MT_ELEMENT_EMPTY,
        MT_COMMENT,
        MT_DECLARATION,
        MT_PROCESSING,
        MT_INVALID
    };

    // Protected methods
    void       init_members(void);
    void       copy_members(const GXml& xml);
    void       free_members(void);
    void       parse(FILE* fptr);
    void       process_markup(GXmlNode** current, const std::string& segment);
    void       process_text(GXmlNode** current, const std::string& segment);
    MarkupType get_markuptype(const std::string& segment) const;

    // Protected data members
    GXmlDocument m_root;         //!< Root node
};

#endif /* GXML_HPP */
