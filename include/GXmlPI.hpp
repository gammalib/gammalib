/***************************************************************************
 *                GXmlPI.hpp - XML PI node class definition                *
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
 * @file GXmlPI.hpp
 * @brief XML PI node class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GXMLPI_HPP
#define GXMLPI_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GUrl.hpp"
#include "GXmlNode.hpp"


/***********************************************************************//**
 * @class GXmlPI
 *
 * @brief XML Processing Instruction node class
 *
 * This class implements a XML Processing Instruction.
 ***************************************************************************/
class GXmlPI : public GXmlNode {

public:
    // Constructors and destructors
    GXmlPI(void);
    GXmlPI(const GXmlPI& node);
    explicit GXmlPI(const std::string& segment);
    virtual ~GXmlPI(void);

    // Operators
    GXmlPI& operator=(const GXmlPI& node);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GXmlPI*     clone(void) const;
    virtual std::string classname(void) const;
    virtual void        write(GUrl& url, const int& indent = 0) const;
    virtual NodeType    type(void) const;
    virtual std::string print(const GChatter& chatter = NORMAL,
                              const int&      indent = 0) const;

    // Other methods
    const std::string&  pi(void) const;
    void                pi(const std::string& pi);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GXmlPI& node);
    void free_members(void);
    void parse(const std::string& segment);

    // Protected data members
    std::string m_pi;       //!< Processing instruction (without brackets)
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GXmlPI").
 ***************************************************************************/
inline
std::string GXmlPI::classname(void) const
{
    return ("GXmlPI");
}


/***********************************************************************//**
 * @brief Return Processing Instruction
 *
 * @return Processing Instruction string.
 ***************************************************************************/
inline
const std::string& GXmlPI::pi(void) const
{
    return (m_pi);
}


/***********************************************************************//**
 * @brief Set Processing Instruction
 *
 * @param[in] pi Processing Instruction string.
 ***************************************************************************/
inline
void GXmlPI::pi(const std::string& pi)
{
    m_pi = pi;
    return;
}


/***********************************************************************//**
 * @brief Return XML node type
 *
 * @return XML node type (NT_PI).
 ***************************************************************************/
inline
GXmlNode::NodeType GXmlPI::type(void) const
{
    return (NT_PI);
}

#endif /* GXMLPI_HPP */
