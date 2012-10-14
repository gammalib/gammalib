/***************************************************************************
 *                GPointing.hpp  -  Pointing interface class               *
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
 * @file GPointing.hpp
 * @brief Pointing interface class definition
 * @author Juergen Knoedlseder
 */

#ifndef GPOINTING_HPP
#define GPOINTING_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GSkyDir.hpp"
#include "GTime.hpp"


/***********************************************************************//**
 * @class GPointing
 *
 * @brief Interface for the pointing classes
 *
 * The pointing class holds information about the time dependent telescope
 * pointing.
 ***************************************************************************/
class GPointing : public GBase {

public:
    // Constructors and destructors
    GPointing(void);
    GPointing(const GPointing& pnt);
    virtual ~GPointing(void);

    // Operators
    virtual GPointing& operator= (const GPointing& pnt);

    // Pure virtual methods
    virtual void           clear(void) = 0;
    virtual GPointing*     clone(void) const = 0;
    virtual const GSkyDir& dir(void) const = 0;
    virtual std::string    print(void) const = 0;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GPointing& pnt);
    void free_members(void);
};

#endif /* GPOINTING_HPP */
