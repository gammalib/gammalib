/***************************************************************************
 *                GXmlPI.hpp - XML PI node class definition                *
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
 * @file GXmlPI.hpp
 * @brief XML PI node class definition
 * @author J. Knodlseder
 */

#ifndef GXMLPI_HPP
#define GXMLPI_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GXmlNode.hpp"


/***********************************************************************//**
 * @class GXmlPI
 *
 * @brief XML PI node class interface defintion.
 *
 * This class implements a XML Processing Instruction.
 ***************************************************************************/
class GXmlPI : public GXmlNode {

public:
    // Constructors and destructors
    GXmlPI(void);
    GXmlPI(const GXmlPI& node);
    GXmlPI(const std::string& segment);
    ~GXmlPI(void);

    // Operators
    GXmlPI& operator= (const GXmlPI& node);

    // Methods
    void     clear(void);
    void     print(std::ostream& os, int indent = 0) const;
    NodeType type(void) const { return NT_PI; }

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GXmlPI& node);
    void         free_members(void);
    GXmlPI*      GXmlPI::clone(void) const;
    void         parse(const std::string& segment);
    std::string  pi(void) const { return m_pi; }
    void         pi(const std::string& pi) { m_pi=pi; }

    // Protected data members
    std::string m_pi;       //!< Processing instruction (without brackets)
};

#endif /* GXMLPI_HPP */
