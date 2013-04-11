/***************************************************************************
 *                       GClass.hpp - [brief descriptor]                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-20xx by [author]                                    *
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
 * @file GClass.hpp
 * @brief [brief descriptor]
 * @author [author]
 */

#ifndef GCLASS_HPP
#define GCLASS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GLog.hpp"


/***********************************************************************//**
 * @class GClass
 *
 * @brief [brief descriptor] interface defintion
 ***************************************************************************/
class GClass {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GClass& c);
    friend GLog&         operator<< (GLog& log, const GClass& c);

public:
    // Constructors and destructors
    GClass(void);
    GClass(const GClass& c);
    virtual ~GClass(void);
 
    // Operators
    GClass& operator= (const GClass& c);

    // Methods
    void        clear(void);
    GClass*     clone(void) const;
    std::string print(void) const;
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GClass& c);
    void free_members(void);

    // Protected data members
    std::string     m_name;          //!< Name
};

#endif /* GCLASS_HPP */
