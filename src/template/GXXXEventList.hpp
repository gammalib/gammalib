/***************************************************************************
 *            GXXXEventList.hpp - [INSTRUMENT] event list class            *
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
 * @file GXXXEventList.hpp
 * @brief [INSTRUMENT] event list class definition
 * @author [AUTHOR]
 */

#ifndef GXXXEVENTLIST_HPP
#define GXXXEVENTLIST_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GEventList.hpp"
#include "GXXXEventAtom.hpp"
#include "GXXXRoi.hpp"

/* __ Forward declarations _______________________________________________ */
class GRoi;
class GFits;
class GFilename;

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GXXXEventList
 *
 * @brief [INSTRUMENT] event list class
 ***************************************************************************/
class GXXXEventList : public GEventList {

public:
    // Constructors and destructors
    GXXXEventList(void);
    explicit GXXXEventList(const GFilename& filename);
    GXXXEventList(const GXXXEventList& list);
    virtual ~GXXXEventList(void);

    // Operators
    virtual GXXXEventList&       operator=(const GXXXEventList& list);
    virtual GXXXEventAtom*       operator[](const int& index);
    virtual const GXXXEventAtom* operator[](const int& index) const;

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GXXXEventList* clone(void) const;
    virtual std::string    classname(void) const;
    virtual int            size(void) const;
    virtual void           load(const GFilename& filename);
    virtual void           save(const GFilename& filename,
                                const bool&      clobber = false) const;
    virtual void           read(const GFits& file);
    virtual void           write(GFits& file) const;
    virtual int            number(void) const;
    virtual void           roi(const GRoi& roi);
    virtual const GXXXRoi& roi(void) const;
    virtual std::string    print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void append(const GXXXEventAtom& event);
    void reserve(const int& number);
    void remove(const int& index, const int& number = 1);

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GXXXEventList& list);
    void         free_members(void);
    virtual void set_energies(void);
    virtual void set_times(void);

    // Protected members
    GXXXRoi                    m_roi;          //!< Region of interest
    std::vector<GXXXEventAtom> m_events;       //!< Events
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GXXXEventList").
 ***************************************************************************/
inline
std::string GXXXEventList::classname(void) const
{
    return ("GXXXEventList");
}


/***********************************************************************//**
 * @brief Return number of events in list
 *
 * @return Number of events in list.
 ***************************************************************************/
inline
int GXXXEventList::size(void) const
{
    return (int)m_events.size();
}


/***********************************************************************//**
 * @brief Return number of events in list
 *
 * @return Number of events in list.
 ***************************************************************************/
inline
int GXXXEventList::number(void) const
{
    return (int)m_events.size();
}


/***********************************************************************//**
 * @brief Return Region of Interest
 *
 * @return Region of Interest.
 ***************************************************************************/
inline
const GXXXRoi& GXXXEventList::roi(void) const
{
    return m_roi;
}


/***********************************************************************//**
 * @brief Reserves space for events
 *
 * @param[in] number Number of events.
 *
 * Reserves space for a @p number events in the event list.
 ***************************************************************************/
inline
void GXXXEventList::reserve(const int& number)
{
    m_events.reserve(number);
    return;
}


/***********************************************************************//**
 * @brief Set energies
 *
 * Sets energies.
 ***************************************************************************/
inline
void GXXXEventList::set_energies(void)
{
    return;
}


/***********************************************************************//**
 * @brief Set times
 *
 * Sets times.
 ***************************************************************************/
inline
void GXXXEventList::set_times(void)
{
    return;
}

#endif /* GLATEVENTLIST_HPP */
