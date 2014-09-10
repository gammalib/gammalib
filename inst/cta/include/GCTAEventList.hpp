/***************************************************************************
 *            GCTAEventList.hpp - CTA event atom container class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
 * @file GCTAEventList.hpp
 * @brief CTA event atom container class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAEVENTLIST_HPP
#define GCTAEVENTLIST_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GEventList.hpp"
#include "GCTAEventAtom.hpp"
#include "GCTARoi.hpp"
#include "GFitsHDU.hpp"
#include "GFitsTable.hpp"
#include "GFitsBinTable.hpp"


/***********************************************************************//**
 * @class GCTAEventList
 *
 * @brief CTA event atom container class
 *
 * This class is a container class for CTA event atoms.
 ***************************************************************************/
class GCTAEventList : public GEventList {

public:
    // Constructors and destructors
    GCTAEventList(void);
    GCTAEventList(const GCTAEventList& list);
    virtual ~GCTAEventList(void);

    // Operators
    virtual GCTAEventList&       operator=(const GCTAEventList& list);
    virtual GCTAEventAtom*       operator[](const int& index);
    virtual const GCTAEventAtom* operator[](const int& index) const;

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GCTAEventList* clone(void) const;
    virtual std::string    classname(void) const;
    virtual int            size(void) const;
    virtual void           load(const std::string& filename);
    virtual void           save(const std::string& filename,
                                const bool& clobber = false) const;
    virtual void           read(const GFits& file);
    virtual void           write(GFits& file) const;
    virtual int            number(void) const;
    virtual void           roi(const GRoi& roi);
    virtual const GCTARoi& roi(void) const;
    std::string            print(const GChatter& chatter = NORMAL) const;

    // Implement other methods
    void   append(const GCTAEventAtom& event);
    void   reserve(const int& number);
    double irf_cache(const std::string& name, const int& index) const;
    void   irf_cache(const std::string& name, const int& index,
                     const double& irf) const;

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GCTAEventList& list);
    void         free_members(void);
    virtual void set_energies(void) { return; }
    virtual void set_times(void) { return; }
    void         read_events(const GFitsTable& hdu);
    void         read_events_v0(const GFitsTable& hdu);
    void         read_events_v1(const GFitsTable& hdu);
    void         read_events_hillas(const GFitsTable& hdu);
    void         read_ds_ebounds(const GFitsHDU& hdu);
    void         read_ds_roi(const GFitsHDU& hdu);
    void         write_events(GFitsBinTable& hdu) const;
    void         write_ds_keys(GFitsHDU& hdu) const;
    int          irf_cache_init(const std::string& name) const;
    int          irf_cache_index(const std::string& name) const;

    // Protected members
    GCTARoi                    m_roi;     //!< Region of interest
    std::vector<GCTAEventAtom> m_events;  //!< Events

    // IRF cache for diffuse models
    mutable std::vector<std::string>          m_irf_names;  //!< Model names
    mutable std::vector<std::vector<double> > m_irf_values; //!< IRF values
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTAEventList").
 ***************************************************************************/
inline
std::string GCTAEventList::classname(void) const
{
    return ("GCTAEventList");
}


/***********************************************************************//**
 * @brief Return number of events in list
 *
 * @return Number of events in list.
 ***************************************************************************/
inline
int GCTAEventList::size(void) const
{
    return (m_events.size());
}


/***********************************************************************//**
 * @brief Return number of events in list
 *
 * @return Number of events in list.
 ***************************************************************************/
inline
int GCTAEventList::number(void) const
{
    return (m_events.size());
}


/***********************************************************************//**
 * @brief Return Region of Interest
 *
 * @return Region of Interest.
 ***************************************************************************/
inline
const GCTARoi& GCTAEventList::roi(void) const
{
    return (m_roi);
}


/***********************************************************************//**
 * @brief Reserves space for events
 *
 * @param[in] number Number of events.
 *
 * Reserves space for number events in the event list.
 ***************************************************************************/
inline
void GCTAEventList::reserve(const int& number)
{
    m_events.reserve(number);
    return;
}

#endif /* GCTAEVENTLIST_HPP */
