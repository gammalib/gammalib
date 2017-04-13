/***************************************************************************
 *                GCTAEventList.hpp - CTA event list class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2017 by Juergen Knoedlseder                         *
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
 * @brief CTA event list class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAEVENTLIST_HPP
#define GCTAEVENTLIST_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GEventList.hpp"
#include "GFilename.hpp"
#include "GCTAEventAtom.hpp"
#include "GCTARoi.hpp"

/* __ Forward declarations _______________________________________________ */
class GFilename;
class GFitsHDU;
class GFitsTable;
class GFitsBinTable;
class GFitsTableCol;

/* __ Constants __________________________________________________________ */
namespace gammalib {
    const std::string extname_cta_events = "EVENTS";
}


/***********************************************************************//**
 * @class GCTAEventList
 *
 * @brief CTA event list class
 *
 * This class implements an event list for CTA.
 ***************************************************************************/
class GCTAEventList : public GEventList {

public:
    // Constructors and destructors
    GCTAEventList(void);
    explicit GCTAEventList(const GFilename& filename);
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
    virtual void           load(const GFilename& filename);
    virtual void           save(const GFilename& filename,
                                const bool&      clobber = false) const;
    virtual void           read(const GFits& fits);
    virtual void           write(GFits& fits) const;
    virtual int            number(void) const;
    virtual void           roi(const GRoi& roi);
    virtual const GCTARoi& roi(void) const;
    std::string            print(const GChatter& chatter = NORMAL) const;

    // Implement other methods
    void               append(const GCTAEventAtom& event);
    void               append_column(const GFitsTableCol& column);
    void               reserve(const int& number);
    void               remove(const int& index, const int& number = 1);
    void               write(GFits& fits, const std::string& evtname,
                                          const std::string& gtiname) const;
    void               fetch(void) const;
    void               dispose(void) const;
    double             irf_cache(const std::string& name, const int& index) const;
    void               irf_cache(const std::string& name, const int& index,
                                 const double& irf) const;
    const std::string& gtiname(void) const;
    void               has_phase(const bool& has_phase);
    void               has_detxy(const bool& has_detxy);
    void               has_mc_id(const bool& has_mc_id);
    const bool&        has_phase() const;
    const bool&        has_detxy() const;
    const bool&        has_mc_id() const;
    
    void               append_column(GFitsTableCol& col);

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GCTAEventList& list);
    void         free_members(void);
    virtual void set_energies(void) { return; }
    virtual void set_times(void) { return; }
    void         read_events(const GFitsTable& table) const;
    void         write_events(GFitsBinTable& table) const;
    void         write_ds_keys(GFitsHDU& hdu,
                               const std::string& gtiname = gammalib::extname_gti) const;
    int          irf_cache_init(const std::string& name) const;
    int          irf_cache_index(const std::string& name) const;

    // Event list meta data
    GCTARoi     m_roi;         //!< Region of interest
    int         m_num_events;  //!< Number of events
    std::string m_gti_extname; //!< GTI extension name
    bool        m_has_phase;   //!< Signal presence of phase
    bool        m_has_detxy;   //!< Signal presence of detector coordinates
    bool        m_has_mc_id;   //!< Signal presence of MC identifier

    // Event list data
    mutable std::vector<GCTAEventAtom>  m_events;   //!< Events
    mutable std::vector<GFitsTableCol*> m_columns;  //!< Pointers to optional columns
    mutable GFilename                   m_filename; //!< Event list file name

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
    return (m_num_events);
}


/***********************************************************************//**
 * @brief Return number of events in list
 *
 * @return Number of events in list.
 ***************************************************************************/
inline
int GCTAEventList::number(void) const
{
    return (m_num_events);
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


/***********************************************************************//**
 * @brief Return Good Time Interval extension name
 *
 * @return Good Time Interval extension name.
 ***************************************************************************/
inline
const std::string& GCTAEventList::gtiname(void) const
{
    return (m_gti_extname);
}


/***********************************************************************//**
 * @brief Set phase information flag
 *
 * @param[in] has_phase Phase information flag.
 *
 * Sets the phase information flag of the event list.
 ***************************************************************************/
inline
void GCTAEventList::has_phase(const bool& has_phase)
{
    m_has_phase = has_phase;
    return;
}


/***********************************************************************//**
 * @brief Set the DETX/DETY information flag
 *
 * @param[in] has_detxy DETX/DETY information flag.
 *
 * Set the DETX/DETY information flag.
 ***************************************************************************/
inline
void GCTAEventList::has_detxy(const bool& has_detxy)
{
    m_has_detxy = has_detxy;
    return;
}


/***********************************************************************//**
 * @brief Set the Monte Carlo identifier flag
 *
 * @param[in] has_mc_id Monte Carlo identifier flag.
 *
 * Set the Monte Carlo identifier flag.
 ***************************************************************************/
inline
void GCTAEventList::has_mc_id(const bool& has_mc_id)
{
    m_has_mc_id = has_mc_id;
    return;
}


/*************************************************************************
 * @brief Return phase information flag
 *
 * @return True if phase information is available, false otherwise.
 ***************************************************************************/
inline
const bool& GCTAEventList::has_phase() const
{
    return m_has_phase;
}


/*************************************************************************
 * @brief Return DETX/DETY information flag
 *
 * @return True if DETX/DETY information is available, false otherwise.
 ***************************************************************************/
inline
const bool& GCTAEventList::has_detxy() const
{
    return m_has_detxy;
}


/*************************************************************************
 * @brief Return Monte Carlo identifier flag
 *
 * @return True if a Monte Carlo identifier is available, false otherwise.
 ***************************************************************************/
inline
const bool& GCTAEventList::has_mc_id() const
{
    return m_has_mc_id;
}

#endif /* GCTAEVENTLIST_HPP */
