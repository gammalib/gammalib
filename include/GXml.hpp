/***************************************************************************
 *                           GXml.hpp - XML class                          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GXml.hpp
 * @brief XML class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GXML_HPP
#define GXML_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GUrl.hpp"
#include "GXmlNode.hpp"
#include "GXmlDocument.hpp"
#include "GXmlElement.hpp"
#include "GXmlText.hpp"

/* __ Forward declarations _______________________________________________ */


/***********************************************************************//**
 * @class GXml
 *
 * @brief XML class
 *
 * This class implements an XML object. It holds the root node of the XML
 * document.
 ***************************************************************************/
class GXml : public GBase {

    // Friend classes
    friend class GXmlNode;
    friend class GXmlDocument;
    friend class GXmlText;

public:
    // Constructors and destructors
    GXml(void);
    GXml(const GXml& xml);
    explicit GXml(const std::string& xml);
    virtual ~GXml(void);

    // Operators
    GXml& operator=(const GXml& xml);

    // Methods
    void         clear(void);
    GXml*        clone(void) const;
    void         append(GXmlNode* node);
    void         load(const std::string& filename);
    void         save(const std::string& filename);
    void         read(GUrl& url);
    void         write(GUrl& url, const int& indent = 0) const;
    int          children(void) const;
    GXmlNode*    child(int index) const;
    int          elements(void) const;
    int          elements(const std::string& name) const;
    GXmlElement* element(int index) const;
    GXmlElement* element(const std::string& name, int index) const;
    std::string  print(void) const;
    std::string  print(const int& indent = 0) const;

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
    void       parse(GUrl& url);
    //int        getchar(GUrl& url, const std::string& string, int& index) const;
    void       process_markup(GXmlNode** current, const std::string& segment);
    void       process_text(GXmlNode** current, const std::string& segment);
    MarkupType get_markuptype(const std::string& segment) const;

    // Protected members
    GXmlDocument m_root;         //!< Root node
};

#endif /* GXML_HPP */
