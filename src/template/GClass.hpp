/***************************************************************************
 *                       GClass.hpp - <brief descriptor>                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-20xx by <author>                                    *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GClass.hpp
 * @brief <brief descriptor>
 * @author <author>
 */

#ifndef GCLASS_HPP
#define GCLASS_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>


/***********************************************************************//**
 * @class GClass
 *
 * @brief <brief descriptor> interface defintion.
 ***************************************************************************/
class GClass {

    // Friend classes
    friend class GModels;
    
    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GClass& c);

public:
    // Constructors and destructors
    GClass(void);
    GClass(const GClass& c);
    ~GClass(void);
 
    // Operators
    GClass& operator= (const GClass& c);

    // Methods
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GClass& c);
    void free_members(void);

    // Protected data members
    std::string     m_name;          //!< Name
};

#endif /* GCLASS_HPP */
