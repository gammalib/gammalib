/***************************************************************************
 *               GEvents.hpp - Abstract event container class              *
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
 * @file GEvents.hpp
 * @brief Abstract event container class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GEVENTS_HPP
#define GEVENTS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GEvent.hpp"
#include "GFits.hpp"
#include "GEbounds.hpp"
#include "GGti.hpp"
#include "GRoi.hpp"


/***********************************************************************//**
 * @class GEvents
 *
 * @brief Abstract event container class
 *
 * This class is an abstract container base class for events. It also holds
 * the energy boundaries and the Good Time Intervals (GTIs) that define the
 * data space. Access to the energy boundaries is provided by the ebounds()
 * methods, the GTI can be accessed by the gti() methods. Furthermore, the
 * emin() and emax() methods return the lowest and highest energy covered
 * by the event container, while the tstart() and tstop() methods return
 * the first and the last time.
 *
 * Access to events (in form of atoms or bins) is provided through the
 * element access operator[].
 *
 * The size() method gives the number of event atoms or bins that is found
 * in the container. The number() method provides the total number of events
 * in the container.
 *
 * An event iterator allows iterating over the event atoms or bins of the
 * container.
 ***************************************************************************/
class GEvents : public GBase {

public:
    // Constructors and destructors
    GEvents(void);
    GEvents(const GEvents& events);
    virtual ~GEvents(void);

    // Operators
    virtual GEvents&      operator=(const GEvents& events);
    virtual GEvent*       operator[](const int& index) = 0;
    virtual const GEvent* operator[](const int& index) const = 0;

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GEvents*    clone(void) const = 0;
    virtual std::string classname(void) const = 0;
    virtual int         size(void) const = 0;
    virtual void        load(const std::string& filename) = 0;
    virtual void        save(const std::string& filename,
                             const bool& clobber = false) const = 0;
    virtual void        read(const GFits& file) = 0;
    virtual void        write(GFits& file) const = 0;
    virtual int         number(void) const = 0;
    virtual std::string print(const GChatter& chatter = NORMAL) const = 0;

    // Implemented methods
    void                ebounds(const GEbounds& ebounds);
    void                gti(const GGti& gti);
    const GEbounds&     ebounds(void) const;
    const GGti&         gti(void) const;
    const GTime&        tstart(void) const;
    const GTime&        tstop(void) const;
    const GEnergy&      emin(void) const;
    const GEnergy&      emax(void) const;

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GEvents& events);
    void         free_members(void);
    virtual void set_energies(void) = 0;
    virtual void set_times(void) = 0;

    // Protected members
    GEbounds m_ebounds;      //!< Energy boundaries covered by events
    GGti     m_gti;          //!< Good time intervals covered by events

};


/***********************************************************************//**
 * @brief Return energy boundaries
 *
 * @return Energy boundaries.
 ***************************************************************************/
inline
const GEbounds& GEvents::ebounds(void) const
{
    return (m_ebounds);
}


/***********************************************************************//**
 * @brief Return Good Time Intervals
 *
 * @return Good Time Intervals.
 ***************************************************************************/
inline
const GGti& GEvents::gti(void) const
{
    return (m_gti);
}


/***********************************************************************//**
 * @brief Return start time
 *
 * @return Start time.
 ***************************************************************************/
inline
const GTime& GEvents::tstart(void) const
{
    return (m_gti.tstart());
}


/***********************************************************************//**
 * @brief Return stop time
 *
 * @return Stop time.
 ***************************************************************************/
inline
const GTime& GEvents::tstop(void) const
{
    return (m_gti.tstop());
}


/***********************************************************************//**
 * @brief Return minimum energy
 *
 * @return Minimum energy.
 ***************************************************************************/
inline
const GEnergy& GEvents::emin(void) const
{
    return (m_ebounds.emin());
}


/***********************************************************************//**
 * @brief Return maximum energy
 *
 * @return Maximum energy.
 ***************************************************************************/
inline
const GEnergy& GEvents::emax(void) const
{
    return (m_ebounds.emax());
}

#endif /* GEVENTS_HPP */
