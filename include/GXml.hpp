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


/***********************************************************************//**
 * @class GXml
 *
 * @brief XML class interface defintion.
 ***************************************************************************/
class GXml {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GXml& xml);

public:
    // Constructors and destructors
    GXml(void);
    GXml(const GXml& xml);
    ~GXml(void);

    // Operators
    GXml& operator= (const GXml& xml);

    // Methods

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GXml& xml);
    void free_members(void);

    // Protected data members
    std::string     m_name;          //!< Name
};

#endif /* GXML_HPP */
