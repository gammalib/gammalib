/***************************************************************************
 *                       GTPLBase.hpp - [WHAT] class                       *
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
 * @file GTPLBase.hpp
 * @brief [WHAT] class definition
 * @author [AUTHOR]
 */

#ifndef GTPLBASE_HPP
#define GTPLBASE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"

/* __ Forward declarations _______________________________________________ */

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GTPLBase
 *
 * @brief [WHAT] class
 *
 * @todo Add class description.
 ***************************************************************************/
class GTPLBase : public GBase {

public:
    // Constructors and destructors
    GTPLBase(void);
    GTPLBase(const GTPLBase& TPL_OBJECT);
    virtual ~GTPLBase(void);

    // Operators
    GTPLBase& operator=(const GTPLBase& TPL_OBJECT);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GTPLBase*   clone(void) const;
    virtual std::string classname(void) const;
    virtual std::string print(const GChatter& chatter = NORMAL) const;

    // Other methods
    // TODO: Add any further methods that are needed

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GTPLBase& TPL_OBJECT);
    void free_members(void);
    
    // Protected members
    // TODO: Add any data members that are necessary
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GTPLBase").
 ***************************************************************************/
inline
std::string GTPLBase::classname(void) const
{
    return ("GTPLBase");
}

#endif /* GTPLBASE_HPP */
