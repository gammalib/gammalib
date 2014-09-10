/***************************************************************************
 *           GXmlDocument.hpp - XML document node class definition         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
 * @file GXmlDocument.hpp
 * @brief XML document node class interface definition
 * @author Juergen Knodlseder
 */

#ifndef GXMLDOCUMENT_HPP
#define GXMLDOCUMENT_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GUrl.hpp"
#include "GXmlNode.hpp"
#include "GXmlAttribute.hpp"


/***********************************************************************//**
 * @class GXmlDocument
 *
 * @brief XML document node class
 *
 * This class implements the root node of an XML document. The root node is
 * a Processing Instruction which contains the following attributes:
 * - @p version
 * - @p encoding
 * - @p standalone
 *
 * All three attributes are systematically written.
 ***************************************************************************/
class GXmlDocument : public GXmlNode {

public:
    // Constructors and destructors
    GXmlDocument(void);
    GXmlDocument(const GXmlDocument& node);
    virtual ~GXmlDocument(void);

    // Operators
    GXmlDocument& operator=(const GXmlDocument& node);

    // Implemented pure virtual base class methods
    virtual void          clear(void);
    virtual GXmlDocument* clone(void) const;
    virtual std::string   classname(void) const;
    virtual void          write(GUrl& url, const int& indent = 0) const;
    virtual NodeType      type(void) const { return NT_DOCUMENT; }
    virtual std::string   print(const GChatter& chatter = NORMAL,
                                const int&      indent = 0) const;

    // Methods
    std::string   version(void) const { return m_version.value(); }
    std::string   encoding(void) const { return m_encoding.value(); }
    std::string   standalone(void) const { return m_standalone.value(); }
    void          version(const std::string& version) { m_version.value(version); }
    void          encoding(const std::string& encoding) { m_encoding.value(encoding); }
    void          standalone(const std::string& standalone) { m_standalone.value(standalone); }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GXmlDocument& node);
    void free_members(void);

    // Protected data members
    GXmlAttribute m_version;      //!< XML version ("1.0", "1.1")
    GXmlAttribute m_encoding;     //!< Encoding (e.g. "UTF-8")
    GXmlAttribute m_standalone;   //!< Standalone ("yes", "no") 
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GXmlDocument").
 ***************************************************************************/
inline
std::string GXmlDocument::classname(void) const
{
    return ("GXmlDocument");
}

#endif /* GXMLDOCUMENT_HPP */
