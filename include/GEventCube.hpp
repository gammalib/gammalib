/***************************************************************************
 *           GEventCube.hpp - Abstract event bin container class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2014 by Juergen Knoedlseder                         *
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
 * @file GEventCube.hpp
 * @brief Abstract event bin container class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GEVENTCUBE_HPP
#define GEVENTCUBE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GEvents.hpp"
#include "GEventBin.hpp"


/***********************************************************************//**
 * @class GEventCube
 *
 * @brief Abstract event bin container class
 *
 * This class is an abstract container class for event bins.
 ***************************************************************************/
class GEventCube : public GEvents {

public:
    // Constructors and destructors
    GEventCube(void);
    GEventCube(const GEventCube& cube);
    virtual ~GEventCube(void);

    // Operators
    virtual GEventCube&      operator=(const GEventCube& cube);
    virtual GEventBin*       operator[](const int& index) = 0;
    virtual const GEventBin* operator[](const int& index) const = 0;

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GEventCube* clone(void) const = 0;
    virtual std::string classname(void) const = 0;
    virtual int         size(void) const = 0;
    virtual int         dim(void) const = 0;
    virtual int         naxis(const int& axis) const = 0;
    virtual void        load(const std::string& filename) = 0;
    virtual void        save(const std::string& filename,
                             const bool& clobber = false) const = 0;
    virtual void        read(const GFits& file) = 0;
    virtual void        write(GFits& file) const = 0;
    virtual int         number(void) const = 0;
    virtual std::string print(const GChatter& chatter = NORMAL) const = 0;

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GEventCube& cube);
    void         free_members(void);
    virtual void set_energies(void);
    virtual void set_times(void);
};


/***********************************************************************//**
 * @brief Set energies (dummy method)
 ***************************************************************************/
inline
void GEventCube::set_energies(void)
{
    return;
}


/***********************************************************************//**
 * @brief Set times (dummy method)
 ***************************************************************************/
inline
void GEventCube::set_times(void)
{
    return;
}

#endif /* GEVENTCUBE_HPP */
