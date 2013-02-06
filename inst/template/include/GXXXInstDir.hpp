/***************************************************************************
 *           GXXXInstDir.hpp  -  XXX instrument direction class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
 * @file GXXXInstDir.hpp
 * @brief XXX instrument direction class definition
 * @author Juergen Knoedlseder
 */

#ifndef GXXXINSTDIR_HPP
#define GXXXINSTDIR_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GInstDir.hpp"


/***********************************************************************//**
 * @class GXXXInstDir
 *
 * @brief Interface for the XXX instrument direction class
 ***************************************************************************/
class GXXXInstDir : public GInstDir {

public:
    // Constructors and destructors
    GXXXInstDir(void);
    GXXXInstDir(const GXXXInstDir& dir);
    virtual ~GXXXInstDir(void);

    // Operators
    GXXXInstDir& operator= (const GXXXInstDir& dir);

    // Methods
    virtual void         clear(void);
    virtual GXXXInstDir* clone(void) const;
    virtual std::string  print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GXXXInstDir& dir);
    void free_members(void);

    // Protected members
};

#endif /* GXXXINSTDIR_HPP */
