/***************************************************************************
 *               GCOMEventList.hpp - COMPTEL event list class              *
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
 * @file GCOMEventList.hpp
 * @brief COMPTEL event list class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCOMEVENTLIST_HPP
#define GCOMEVENTLIST_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GEventList.hpp"
#include "GCOMEventAtom.hpp"
#include "GCOMRoi.hpp"

/* __ Forward declarations _______________________________________________ */
class GRoi;
class GFits;
class GFilename;

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GCOMEventList
 *
 * @brief COMPTEL event list class
 ***************************************************************************/
class GCOMEventList : public GEventList {

public:
    // Constructors and destructors
    GCOMEventList(void);
    explicit GCOMEventList(const GFilename& filename);
    GCOMEventList(const GCOMEventList& list);
    virtual ~GCOMEventList(void);

    // Operators
    virtual GCOMEventList&       operator=(const GCOMEventList& list);
    virtual GCOMEventAtom*       operator[](const int& index);
    virtual const GCOMEventAtom* operator[](const int& index) const;

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GCOMEventList* clone(void) const;
    virtual std::string    classname(void) const;
    virtual int            size(void) const;
    virtual void           load(const GFilename& filename);
    virtual void           save(const GFilename& filename,
                                const bool&      clobber = false) const;
    virtual void           read(const GFits& file);
    virtual void           write(GFits& file) const;
    virtual int            number(void) const;
    virtual void           roi(const GRoi& roi);
    virtual const GCOMRoi& roi(void) const;
    virtual std::string    print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void append(const GCOMEventAtom& event);
    void reserve(const int& number);
    void remove(const int& index, const int& number = 1);

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GCOMEventList& list);
    void         free_members(void);
    virtual void set_energies(void);
    virtual void set_times(void);

    // Protected members
    GCOMRoi                    m_roi;          //!< Region of interest
    std::vector<GCOMEventAtom> m_events;       //!< Events
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMEventList").
 ***************************************************************************/
inline
std::string GCOMEventList::classname(void) const
{
    return ("GCOMEventList");
}


/***********************************************************************//**
 * @brief Return number of events in list
 *
 * @return Number of events in list.
 ***************************************************************************/
inline
int GCOMEventList::size(void) const
{
    return (int)m_events.size();
}


/***********************************************************************//**
 * @brief Return number of events in list
 *
 * @return Number of events in list.
 ***************************************************************************/
inline
int GCOMEventList::number(void) const
{
    return (int)m_events.size();
}


/***********************************************************************//**
 * @brief Return Region of Interest
 *
 * @return Region of Interest.
 ***************************************************************************/
inline
const GCOMRoi& GCOMEventList::roi(void) const
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
void GCOMEventList::reserve(const int& number)
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
void GCOMEventList::set_energies(void)
{
    return;
}


/***********************************************************************//**
 * @brief Set times
 *
 * Sets times.
 ***************************************************************************/
inline
void GCOMEventList::set_times(void)
{
    return;
}

#endif /* GLATEVENTLIST_HPP */
