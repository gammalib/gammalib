/***************************************************************************
 *                  GCOMDri.hpp - COMPTEL Data Space class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knoedlseder                              *
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
 * @file GCOMDri.hpp
 * @brief COMPTEL Data Space class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCOMDRI_HPP
#define GCOMDRI_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"

/* __ Forward declarations _______________________________________________ */

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GCOMDri
 *
 * @brief COMPTEL Data Space class
 *
 * @todo Add class description.
 ***************************************************************************/
class GCOMDri : public GBase {

public:
    // Constructors and destructors
    GCOMDri(void);
    GCOMDri(const GCOMDri& dri);
    virtual ~GCOMDri(void);

    // Operators
    GCOMDri& operator=(const GCOMDri& dri);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GCOMDri*   clone(void) const;
    virtual std::string classname(void) const;
    virtual std::string print(const GChatter& chatter = NORMAL) const;

    // Other methods
    // TODO: Add any further methods that are needed

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMDri& dri);
    void free_members(void);
    
    // Protected members
    // TODO: Add any data members that are necessary
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMDri").
 ***************************************************************************/
inline
std::string GCOMDri::classname(void) const
{
    return ("GCOMDri");
}

#endif /* GCOMDRI_HPP */
