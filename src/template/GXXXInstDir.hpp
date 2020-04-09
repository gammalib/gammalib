/***************************************************************************
 *         GXXXInstDir.hpp - [INSTRUMENT] instrument direction class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) [YEAR] by [AUTHOR]                                       *
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
 * @brief [INSTRUMENT] instrument direction class definition
 * @author [AUTHOR]
 */

#ifndef GXXXINSTDIR_HPP
#define GXXXINSTDIR_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GInstDir.hpp"

/* __ Forward declarations _______________________________________________ */

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GXXXInstDir
 *
 * @brief [INSTRUMENT] instrument direction class
 *
 * The [INSTRUMENT] instrument direction defines the spatial information
 * associated to an event.
 ***************************************************************************/
class GXXXInstDir : public GInstDir {

public:
    // Constructors and destructors
    GXXXInstDir(void);
    GXXXInstDir(const GXXXInstDir& dir);
    virtual ~GXXXInstDir(void);

    // Operators
    GXXXInstDir& operator=(const GXXXInstDir& dir);

    // Implemented pure virtual base class methods
    virtual void         clear(void);
    virtual GXXXInstDir* clone(void) const;
    virtual std::string  classname(void) const;
    virtual double       hash(void) const;
    virtual std::string  print(const GChatter& chatter = NORMAL) const;

    // Other methods
    // TODO: Add any further methods that are needed

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GXXXInstDir& dir);
    void free_members(void);
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GXXXInstDir").
 ***************************************************************************/
inline
std::string GXXXInstDir::classname(void) const
{
    return ("GXXXInstDir");
}


/***********************************************************************//**
 * @brief Return instrument direction hash value
 *
 * @return Hash value.
 *
 * Returns a hash value that can be used in the response cache.
 ***************************************************************************/
inline
double GInstDir::hash(void) const
{
    // TODO: Remplace by a reasonable hash algorithm, assuring that values
    // are comprised between 0 and 1e6
    return 12345.67;
}

#endif /* GXXXINSTDIR_HPP */
