/***************************************************************************
 *                GXmlNode.hpp - Abstract XML node base class              *
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
 * @file GXmlNode.hpp
 * @brief Abstract XML node base class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GXMLNODE_HPP
#define GXMLNODE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GBase.hpp"


/***********************************************************************//**
 * @class GXmlNode
 *
 * @brief Abstract XML node base class
 *
 * This class defines an abstract node of a XML document.
 ***************************************************************************/
class GXmlNode : public GBase {

public:
    // Constructors and destructors
    GXmlNode(void);
    GXmlNode(const GXmlNode& node);
    virtual ~GXmlNode(void);

    // Operators
    GXmlNode& operator= (const GXmlNode& node);

    // Public enumerators
    enum NodeType {
        NT_DOCUMENT,
        NT_ELEMENT,
        NT_COMMENT,
        NT_UNKNOWN,
        NT_TEXT,
        NT_DECLARATION,
        NT_PI,
        NT_TYPECOUNT
    };

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GXmlNode*   clone(void) const = 0;
    virtual void        write(FILE* fptr, int indent = 0) const = 0;
    virtual NodeType    type(void) const = 0;
    virtual std::string print(void) const;
    virtual std::string print(int indent) const = 0;

    // Methods
    void      append(GXmlNode* node);
    int       children(void) const;
    GXmlNode* child(int index) const;
    int       elements(void) const;
    int       elements(const std::string& name) const;
    GXmlNode* element(int index) const;
    GXmlNode* element(const std::string& name, int index) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GXmlNode& node);
    void free_members(void);

    // Protected data members
    std::vector<GXmlNode*> m_nodes;    //!< Pointer to nodes contained in node
};

#endif /* GXMLNODE_HPP */
