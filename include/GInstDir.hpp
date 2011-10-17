/***************************************************************************
 *         GInstDir.hpp  -  Instrument direction abstract base class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
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
 * @file GInstDir.hpp
 * @brief GInstDir abstract base class definition.
 * @author J. Knodlseder
 */

#ifndef GINSTDIR_HPP
#define GINSTDIR_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GLog.hpp"


/***********************************************************************//**
 * @class GInstDir
 *
 * @brief Abstract interface for the instrument direction of an event.
 *
 * The instrument direction of an event is the equivalent of the sky
 * direction (implemented by GSkyDir) but in the instrument data space.
 * The instrument direction may be any kind of position or direction
 * information encoded in the data space, such as incident event
 * reconstructions for imaging devices or detector numbers etc. for
 * non-imaging devices.
 ***************************************************************************/
class GInstDir {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GInstDir& dir);
    friend GLog&         operator<< (GLog& log, const GInstDir& dir);

public:
    // Constructors and destructors
    GInstDir(void);
    GInstDir(const GInstDir& dir);
    virtual ~GInstDir(void);

    // Operators
    virtual GInstDir& operator= (const GInstDir& dir);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GInstDir*   clone(void) const = 0;
    virtual std::string print(void) const = 0;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GInstDir& dir);
    void free_members(void);
};

#endif /* GINSTDIR_HPP */
