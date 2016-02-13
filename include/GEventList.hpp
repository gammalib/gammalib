/***************************************************************************
 *            GEventList.hpp - Abstract event atom container class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2016 by Juergen Knoedlseder                         *
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
 * @file GEventList.hpp
 * @brief Abstract event atom container class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GEVENTLIST_HPP
#define GEVENTLIST_HPP

/* __ Includes ___________________________________________________________ */
#include "GEvents.hpp"
#include "GEventAtom.hpp"

/* __ Forward declarations _______________________________________________ */
class GFilename;


/***********************************************************************//**
 * @class GEventList
 *
 * @brief Abstract event atom container class
 *
 * This class is an abstract container class for event atoms.
 *
 * In addition to the base class methods, it defines also an interface to the
 * region of interest (ROI) that is covered by the data space. The ROI can
 * be accessed by the roi() methods.
 ***************************************************************************/
class GEventList : public GEvents {

public:
    // Constructors and destructors
    GEventList(void);
    GEventList(const GEventList& list);
    virtual ~GEventList(void);

    // Operators
    virtual GEventList&       operator=(const GEventList& list);
    virtual GEventAtom*       operator[](const int& index) = 0;
    virtual const GEventAtom* operator[](const int& index) const = 0;

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GEventList* clone(void) const = 0;
    virtual std::string classname(void) const = 0;
    virtual int         size(void) const = 0;
    virtual void        load(const GFilename& filename) = 0;
    virtual void        save(const GFilename& filename,
                             const bool&      clobber = false) const = 0;
    virtual void        read(const GFits& file) = 0;
    virtual void        write(GFits& file) const = 0;
    virtual int         number(void) const = 0;
    virtual void        roi(const GRoi& roi) = 0;
    virtual const GRoi& roi(void) const = 0;
    virtual std::string print(const GChatter& chatter = NORMAL) const = 0;

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GEventList& list);
    void         free_members(void);
    virtual void set_energies(void) = 0;
    virtual void set_times(void) = 0;
};

#endif /* GEVENTLIST_HPP */
