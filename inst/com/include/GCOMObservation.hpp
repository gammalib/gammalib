/***************************************************************************
 *            GCOMObservation.hpp - COMPTEL observation class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2023 by Juergen Knoedlseder                         *
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
 * @file GCOMObservation.hpp
 * @brief COMPTEL observation class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCOMOBSERVATION_HPP
#define GCOMOBSERVATION_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GObservation.hpp"
#include "GTime.hpp"
#include "GFilename.hpp"
#include "GCOMResponse.hpp"
#include "GCOMTim.hpp"
#include "GCOMOads.hpp"
#include "GCOMHkds.hpp"
#include "GCOMBvcs.hpp"
#include "GCOMDri.hpp"
#include "GCOMEventList.hpp"
#include "GCOMEventCube.hpp"

/* __ Forward declarations _______________________________________________ */
class GCaldb;
class GResponse;
class GModelSky;
class GModels;
class GXmlElement;
class GFitsHDU;
class GSkyMap;
class GCOMStatus;


/***********************************************************************//**
 * @class GCOMObservation
 *
 * @brief Interface class for COMPTEL observations
 *
 * This class implements a COMPTEL observation. Each COMPTEL observation is
 * defined for a given energy range, and is composed of a DRE, DRB, DRG and
 * DRX file. The DRE file contains the event data, the DRB file contains a
 * background model, the DRG file contains geometry factors, and the DRX file
 * contains the exposure.
 ***************************************************************************/
class GCOMObservation : public GObservation {

public:
    // Constructors and destructors
    GCOMObservation(void);
    explicit GCOMObservation(const GXmlElement& xml);
    GCOMObservation(const GCOMDri& dre,
                    const GCOMDri& drb,
                    const GCOMDri& drg,
                    const GCOMDri& drx);
    GCOMObservation(const GCOMDri& dre,
                    const GCOMDri& drb,
                    const GCOMDri& drw,
                    const GCOMDri& drg,
                    const GCOMDri& drx);
    GCOMObservation(const GFilename& drename,
                    const GFilename& drbname,
                    const GFilename& drwname,
                    const GFilename& drgname,
                    const GFilename& drxname);
    GCOMObservation(const GFilename&              evpname,
                    const GFilename&              timname,
                    const std::vector<GFilename>& oadnames,
                    const std::vector<GFilename>& hkdnames = std::vector<GFilename> (),
                    const GFilename&              bvcname = "");
    GCOMObservation(const GCOMObservation& obs);
    virtual ~GCOMObservation(void);

    // Operators
    virtual GCOMObservation& operator=(const GCOMObservation& obs);

    // Implement pure virtual methods
    virtual void                clear(void);
    virtual GCOMObservation*    clone(void) const;
    virtual std::string         classname(void) const;
    virtual void                response(const GResponse& rsp);
    virtual const GCOMResponse* response(void) const;
    virtual std::string         instrument(void) const;
    virtual double              ontime(void) const;
    virtual double              livetime(void) const;
    virtual double              deadc(const GTime& time = GTime()) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;
    virtual std::string         print(const GChatter& chatter = NORMAL) const;

    // Reimplemented virtual methods
    virtual double   npred(const GModel& model) const;

    // Other methods
    bool             is_unbinned(void) const;
    bool             is_binned(void) const;
    void             load(const GFilename& drename,
                          const GFilename& drbname,
                          const GFilename& drwname,
                          const GFilename& drgname,
                          const GFilename& drxname);
    void             load(const GFilename&              evpname,
                          const GFilename&              timname,
                          const std::vector<GFilename>& oadnames,
                          const std::vector<GFilename>& hkdnames = std::vector<GFilename> (),
                          const GFilename&              bvcname = "");
    void             response(const GCaldb& caldb, const std::string& rspname);
    void             response(const GCOMResponse& response);
    void             obs_id(const double& id);
    void             ontime(const double& ontime);
    void             livetime(const double& livetime);
    void             deadc(const double& deadc);
    void             ewidth(const double& ewidth);
    const double&    obs_id(void) const;
    const double&    ewidth(void) const;
    const GCOMDri&   drb(void) const;
    const GCOMDri&   drw(void) const;
    const GCOMDri&   drg(void) const;
    const GCOMDri&   drx(void) const;
    GCOMDri          drm(const GModels& models) const;
    const GCOMTim&   tim(void) const;
    void             tim(const GCOMTim& tim);
    const GCOMOads&  oads(void) const;
    void             oads(const GCOMOads& oads);
    const GCOMHkds&  hkds(void) const;
    void             hkds(const GCOMHkds& hkds);
    const GCOMBvcs&  bvcs(void) const;
    void             bvcs(const GCOMBvcs& bvcs);
    const GFilename& drename(void) const;
    const GFilename& drbname(void) const;
    const GFilename& drwname(void) const;
    const GFilename& drgname(void) const;
    const GFilename& drxname(void) const;
    const GFilename& rspname(void) const;
    const int&       phi_first(void) const;
    const int&       phi_last(void) const;
    void             drename(const GFilename& drename);
    void             drbname(const GFilename& drbname);
    void             drwname(const GFilename& drwname);
    void             drgname(const GFilename& drgname);
    void             drxname(const GFilename& drxname);
    void             rspname(const GFilename& rspname);
    void             phi_first(const int& phi_first);
    void             phi_last(const int& phi_last);
    void             compute_drb(const std::string& method,
                                 const GCOMDri&     drm,
                                 const int&         nrunav = 3,
                                 const int&         navgr  = 3,
                                 const int&         nincl  = 15,
                                 const int&         nexcl  = 0);

protected:
    // Protected methods
    void    init_members(void);
    void    copy_members(const GCOMObservation& obs);
    void    free_members(void);
    void    load_dre(const GFilename& drename);
    void    load_drb(const GFilename& drbname);
    void    load_drw(const GFilename& drwname);
    void    load_drg(const GFilename& drgname);
    void    load_drx(const GFilename& drxname);
    bool    check_dri(const GCOMDri& map) const;
    void    read_attributes(const GFitsHDU* hdu);
    void    write_attributes(GFitsHDU* hdu) const;
    void    compute_drb_phinor(const GCOMDri& drm);
    void    compute_drb_bgdlixa(const GCOMDri& drm,
                                const int&     nrunav = 3,
                                const int&     navgr  = 3,
                                const int&     nincl  = 15,
                                const int&     nexcl  = 0);
    void    compute_drb_bgdlixe(const GCOMDri& drm,
                                const int&     nrunav = 3,
                                const int&     navgr  = 3,
                                const int&     nincl  = 15,
                                const int&     nexcl  = 0);
    void    compute_drb_bgdlixf(const GCOMDri& drm,
                                const int&     nrunav = 3,
                                const int&     navgr  = 3,
                                const int&     nincl  = 15,
                                const int&     nexcl  = 0);
    GSkyMap get_weighted_drg_map(void) const;
    void    get_bgdlixa_phibar_indices(const int& iphibar,
                                       const int& nincl,
                                       const int& nexcl,
                                       int*       isel1,
                                       int*       iex1,
                                       int*       iex2,
                                       int*       isel2) const;

    // Overwritten virtual methods
    virtual bool use_event_for_likelihood(const int& index) const;

    // Protected members
    std::string            m_instrument; //!< Instrument name
    GCOMResponse           m_response;   //!< Response functions
    double                 m_obs_id;     //!< Observation ID
    double                 m_ontime;     //!< Ontime (sec)
    double                 m_livetime;   //!< Livetime (sec)
    double                 m_deadc;      //!< Deadtime correction

    // Protected members for binned observation
    GFilename              m_drename;    //!< DRE filename
    GFilename              m_drbname;    //!< DRB filename
    GFilename              m_drwname;    //!< DRW filename
    GFilename              m_drgname;    //!< DRG filename
    GFilename              m_drxname;    //!< DRX filename
    GFilename              m_rspname;    //!< Response cache filename
    GCOMDri                m_drb;        //!< Background model
    GCOMDri                m_drw;        //!< Weighting cube
    GCOMDri                m_drg;        //!< Geometry factors
    GCOMDri                m_drx;        //!< Exposure map
    double                 m_ewidth;     //!< Energy width (MeV)
    int                    m_phi_first;  //!< First Phibar layer to use for likelihood
    int                    m_phi_last;   //!< Last Phibar layer to use for likelihood

    // Protected members for unbinned observation
    GFilename              m_evpname;    //!< EVP filename
    GFilename              m_timname;    //!< TIM filename
    std::vector<GFilename> m_oadnames;   //!< OAD filenames
    std::vector<GFilename> m_hkdnames;   //!< HKD filenames
    GFilename              m_bvcname;    //!< BVC filename
    GCOMTim                m_tim;        //!< COMPTEL Good Time Intervals
    GCOMOads               m_oads;       //!< Orbit Aspect Data
    GCOMHkds               m_hkds;       //!< Housekeeping Data
    GCOMBvcs               m_bvcs;       //!< Solar System Barycentre Data
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMObservation").
 ***************************************************************************/
inline
std::string GCOMObservation::classname(void) const
{
    return ("GCOMObservation");
}


/***********************************************************************//**
 * @brief Return response function
 *
 * @return Response function.
 ***************************************************************************/
inline
const GCOMResponse* GCOMObservation::response(void) const
{
    // Return response pointer
    return &m_response;
}


/***********************************************************************//**
 * @brief Set response function
 *
 * @param[in] response Response function.
 ***************************************************************************/
inline
void GCOMObservation::response(const GCOMResponse& response)
{
    m_response = response;
    return;
}


/***********************************************************************//**
 * @brief Return instrument
 *
 * @return Instrument name.
 ***************************************************************************/
inline
std::string GCOMObservation::instrument(void) const
{
    // Return instrument
    return (m_instrument);
}


/***********************************************************************//**
 * @brief Return ontime
 *
 * @return Ontime (seconds).
 ***************************************************************************/
inline
double GCOMObservation::ontime(void) const
{
    // Return ontime
    return (m_ontime);
}


/***********************************************************************//**
 * @brief Return livetime
 *
 * @return Livetime (seconds).
 ***************************************************************************/
inline
double GCOMObservation::livetime(void) const
{
    // Return livetime
    return (m_livetime);
}


/***********************************************************************//**
 * @brief Return deadtime correction factor
 *
 * @param[in] time Time.
 *
 * @return Deadtime correction factor.
 ***************************************************************************/
inline
double GCOMObservation::deadc(const GTime& time) const
{
    // Return livetime
    return (m_deadc);
}


/***********************************************************************//**
 * @brief Set observation ID
 *
 * @param[in] id Observation ID.
 ***************************************************************************/
inline
void GCOMObservation::obs_id(const double& id)
{
    m_obs_id = id;
    return;
}


/***********************************************************************//**
 * @brief Set ontime
 *
 * @param[in] ontime Ontime.
 ***************************************************************************/
inline
void GCOMObservation::ontime(const double& ontime)
{
    m_ontime = ontime;
    return;
}


/***********************************************************************//**
 * @brief Set livetime
 *
 * @param[in] livetime Livetime.
 ***************************************************************************/
inline
void GCOMObservation::livetime(const double& livetime)
{
    m_livetime = livetime;
    return;
}


/***********************************************************************//**
 * @brief Set deadtime correction factor
 *
 * @param[in] deadc Deadtime correction factor.
 ***************************************************************************/
inline
void GCOMObservation::deadc(const double& deadc)
{
    m_deadc = deadc;
    return;
}


/***********************************************************************//**
 * @brief Set energy width
 *
 * @param[in] ewidth Energy width (MeV).
 ***************************************************************************/
inline
void GCOMObservation::ewidth(const double& ewidth)
{
    m_ewidth = ewidth;
    return;
}


/***********************************************************************//**
 * @brief Return observation ID
 *
 * @return Observation ID.
 ***************************************************************************/
inline
const double& GCOMObservation::obs_id(void) const
{
    // Return observation ID
    return (m_obs_id);
}


/***********************************************************************//**
 * @brief Return energy width
 *
 * @return Energy width (MeV).
 ***************************************************************************/
inline
const double& GCOMObservation::ewidth(void) const
{
    // Return energy width
    return (m_ewidth);
}


/***********************************************************************//**
 * @brief Return background model
 *
 * @return Background model.
 ***************************************************************************/
inline
const GCOMDri& GCOMObservation::drb(void) const
{
    // Return background model
    return (m_drb);
}


/***********************************************************************//**
 * @brief Return weighting cube
 *
 * @return Weighting cube.
 ***************************************************************************/
inline
const GCOMDri& GCOMObservation::drw(void) const
{
    // Return weighting cube
    return (m_drw);
}


/***********************************************************************//**
 * @brief Return geometry factors
 *
 * @return Geometry factors.
 ***************************************************************************/
inline
const GCOMDri& GCOMObservation::drg(void) const
{
    // Return geometry factors
    return (m_drg);
}


/***********************************************************************//**
 * @brief Return exposure
 *
 * @return Exposure.
 ***************************************************************************/
inline
const GCOMDri& GCOMObservation::drx(void) const
{
    // Return exposure
    return (m_drx);
}


/***********************************************************************//**
 * @brief Return COMPTEL Good Time Intervals
 *
 * @return COMPTEL Good Time Intervals.
 ***************************************************************************/
inline
const GCOMTim& GCOMObservation::tim(void) const
{
    // Return COMPTEL Good Time Intervals
    return (m_tim);
}


/***********************************************************************//**
 * @brief Set COMPTEL Good Time Intervals
 *
 * @param[in] tim COMPTEL Good Time Intervals.
 ***************************************************************************/
inline
void GCOMObservation::tim(const GCOMTim& tim)
{
    m_tim = tim;
    return;
}


/***********************************************************************//**
 * @brief Return Orbit Aspect Data
 *
 * @return Orbit Aspect Data
 ***************************************************************************/
inline
const GCOMOads& GCOMObservation::oads(void) const
{
    // Return Orbit Aspect Data
    return (m_oads);
}


/***********************************************************************//**
 * @brief Set Orbit Aspect Data
 *
 * @param[in] oads Orbit Aspect Data.
 ***************************************************************************/
inline
void GCOMObservation::oads(const GCOMOads& oads)
{
    m_oads = oads;
    return;
}


/***********************************************************************//**
 * @brief Return Housekeeping Data collection
 *
 * @return Housekeeping Data collection
 ***************************************************************************/
inline
const GCOMHkds& GCOMObservation::hkds(void) const
{
    // Return Housekeeping Data collection
    return (m_hkds);
}


/***********************************************************************//**
 * @brief Set Housekeeping Data collection
 *
 * @param[in] hkds Housekeeping Data collection.
 ***************************************************************************/
inline
void GCOMObservation::hkds(const GCOMHkds& hkds)
{
    m_hkds = hkds;
    return;
}


/***********************************************************************//**
 * @brief Return Solar System Barycentre Data
 *
 * @return Solar System Barycentre Data
 ***************************************************************************/
inline
const GCOMBvcs& GCOMObservation::bvcs(void) const
{
    // Return Solar System Barycentre Data
    return (m_bvcs);
}


/***********************************************************************//**
 * @brief Set Solar System Barycentre Data
 *
 * @param[in] bvcs Solar System Barycentre Data.
 ***************************************************************************/
inline
void GCOMObservation::bvcs(const GCOMBvcs& bvcs)
{
    m_bvcs = bvcs;
    return;
}


/***********************************************************************//**
 * @brief Check whether observation is unbinned
 *
 * @return True if observation is unbinned.
 ***************************************************************************/
inline
bool GCOMObservation::is_unbinned(void) const
{
    return (dynamic_cast<const GCOMEventList*>(m_events) != NULL);
}


/***********************************************************************//**
 * @brief Check whether observation is binned
 *
 * @return True if observation is unbinned.
 ***************************************************************************/
inline
bool GCOMObservation::is_binned(void) const
{
    return (dynamic_cast<const GCOMEventCube*>(m_events) != NULL);
}


/***********************************************************************//**
 * @brief Return DRE filename
 *
 * @return DRE filename.
 ***************************************************************************/
inline
const GFilename& GCOMObservation::drename(void) const
{
    // Return DRE filename
    return (m_drename);
}


/***********************************************************************//**
 * @brief Return DRB filename
 *
 * @return DRB filename.
 ***************************************************************************/
inline
const GFilename& GCOMObservation::drbname(void) const
{
    // Return DRB filename
    return (m_drbname);
}


/***********************************************************************//**
 * @brief Return DRW filename
 *
 * @return DRW filename.
 ***************************************************************************/
inline
const GFilename& GCOMObservation::drwname(void) const
{
    // Return DRW filename
    return (m_drwname);
}


/***********************************************************************//**
 * @brief Return DRG filename
 *
 * @return DRG filename.
 ***************************************************************************/
inline
const GFilename& GCOMObservation::drgname(void) const
{
    // Return DRG filename
    return (m_drgname);
}


/***********************************************************************//**
 * @brief Return DRX filename
 *
 * @return DRX filename.
 ***************************************************************************/
inline
const GFilename& GCOMObservation::drxname(void) const
{
    // Return DRX filename
    return (m_drxname);
}


/***********************************************************************//**
 * @brief Return response cache filename
 *
 * @return Response cache filename.
 ***************************************************************************/
inline
const GFilename& GCOMObservation::rspname(void) const
{
    // Return response cache filename
    return (m_rspname);
}


/***********************************************************************//**
 * @brief Return index of first Phibar layer to be used for likelihood fitting
 *
 * @return Index of first Phibar layer.
 ***************************************************************************/
inline
const int& GCOMObservation::phi_first(void) const
{
    // Return index of first Phibar layer
    return (m_phi_first);
}


/***********************************************************************//**
 * @brief Return index of last Phibar layer to be used for likelihood fitting
 *
 * @return Index of last Phibar layer.
 ***************************************************************************/
inline
const int& GCOMObservation::phi_last(void) const
{
    // Return index of last Phibar layer
    return (m_phi_last);
}


/***********************************************************************//**
 * @brief Set DRE filename
 *
 * @param[in] drename DRE filename.
 ***************************************************************************/
inline
void GCOMObservation::drename(const GFilename& drename)
{
    m_drename = drename;
    return;
}


/***********************************************************************//**
 * @brief Set DRB filename
 *
 * @param[in] drbname DRB filename.
 ***************************************************************************/
inline
void GCOMObservation::drbname(const GFilename& drbname)
{
    m_drbname = drbname;
    return;
}


/***********************************************************************//**
 * @brief Set DRW filename
 *
 * @param[in] drwname DRW filename.
 ***************************************************************************/
inline
void GCOMObservation::drwname(const GFilename& drwname)
{
    m_drwname = drwname;
    return;
}


/***********************************************************************//**
 * @brief Set DRG filename
 *
 * @param[in] drgname DRG filename.
 ***************************************************************************/
inline
void GCOMObservation::drgname(const GFilename& drgname)
{
    m_drgname = drgname;
    return;
}


/***********************************************************************//**
 * @brief Set DRX filename
 *
 * @param[in] drxname DRX filename.
 ***************************************************************************/
inline
void GCOMObservation::drxname(const GFilename& drxname)
{
    m_drxname = drxname;
    return;
}


/***********************************************************************//**
 * @brief Set response cache filename filename
 *
 * @param[in] cachename Response cache filename.
 ***************************************************************************/
inline
void GCOMObservation::rspname(const GFilename& rspname)
{
    m_rspname = rspname;
    return;
}


/***********************************************************************//**
 * @brief Set index of first Phibar layer to be used for likelihood fitting
 *
 * @param[in] phi_first Index of first Phibar layer.
 ***************************************************************************/
inline
void GCOMObservation::phi_first(const int& phi_first)
{
    m_phi_first = phi_first;
    return;
}


/***********************************************************************//**
 * @brief Set index of last Phibar layer to be used for likelihood fitting
 *
 * @param[in] phi_last Index of last Phibar layer.
 ***************************************************************************/
inline
void GCOMObservation::phi_last(const int& phi_last)
{
    m_phi_last = phi_last;
    return;
}

#endif /* GCOMOBSERVATION_HPP */
